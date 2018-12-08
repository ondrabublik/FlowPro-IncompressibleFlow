package flowpro.user.equation;

import flowpro.api.ElementData;
import flowpro.api.Equation;
import flowpro.api.FlowProProperties;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

public class IncompressibleNavierStokes implements Equation {

    protected class BoundaryType {

        static final int WALL = -1;
        static final int INLET = -2;
        static final int OUTLET = -3;
        static final int INVISCID_WALL = -4;
    }

    protected int dim;
    protected int nEqs;
    protected boolean isDiffusive;

    protected double Re; // Reynoldsovo cislo

    // reference values
    protected double lRef;
    protected double pRef;
    protected double rhoRef;
    protected double velocityRef;
    protected double tRef;

    // inlet boundary condition
    protected boolean isInletSupersonic;
    // subsonic inlet boundary condition 
    protected final double pIn0 = 1; // static pressure
    protected final double rhoIn0 = 1; // static density
    // supersonic inlet boundary condition
    protected double[] VIn;

    // outlet boundary condition
    protected double pOut; // pressure

    double gravityAcceleration;

    @Override
    public int dim() {
        return dim;
    }

    @Override
    public int nEqs() {
        return nEqs;
    }

    @Override
    public boolean isConvective() {
        return true;
    }

    @Override
    public boolean isDiffusive() {
        return true;
    }

    @Override
    public void init(FlowProProperties props) throws IOException {

        dim = props.getInt("dimension");
        nEqs = dim + 1;

        gravityAcceleration = 9.81;

        isDiffusive = props.getBoolean("isFlowViscous");

        // reference values from inlet
        VIn = props.getDoubleArray("VIn");

        // other reference values
        if (props.containsKey("lRef")) {
            lRef = props.getDouble("lRef");
        } else {
            lRef = 1;
        }

        velocityRef = 0;
        for (int d = 0; d < dim; ++d) {
            velocityRef += VIn[d] * VIn[d];
        }
        velocityRef = Math.sqrt(velocityRef);
        tRef = lRef / velocityRef;

        // outlet
        if (props.containsKey("pOut")) {
            pOut = props.getDouble("pOut");
        } else {
            throw new IOException("outlet boundary pressure condition is not specified");
        }

        Re = props.getDouble("reynolds");
    }

    @Override
    public void setState(double dt, double t) {
    }

    @Override
    public double[] constInitCondition() {
        if (dim == 2) {
            return new double[]{pOut, VIn[0], VIn[1]};
        } else {
            return new double[]{pOut, VIn[0], VIn[1], VIn[2]};
        }
    }

    @Override
    public void limitUnphysicalValues(double[] Ws, double[] W, int nBasis) { // limituje zaporne hodnoty
    }

    //  nevazky tok stenou _____________________________________________________
    @Override
    public double[] numericalConvectiveFlux(double WL[], double WR[], double[] n, int TT, ElementData elem) {
        double[] f = new double[nEqs];

        switch (TT) {
            case (BoundaryType.WALL):
            case (BoundaryType.INVISCID_WALL):
                f[0] = 0;
                for (int d = 0; d < dim; ++d) {
                    f[d + 1] = WL[0] * n[d];
                }
                break;
            case (BoundaryType.INLET):
            case (BoundaryType.OUTLET):
                f = convectiveFlux(WR, n, elem);
                break;

            default: // vnitrni stena
//                double[] fL = convectiveFlux(WL, n, elem);
//                double[] fR = convectiveFlux(WR, n, elem);
//                double maxEigenValue = Math.max(maxEigenvalue(WL, elem), maxEigenvalue(WR, elem));
//                for (int j = 0; j < nEqs; j++) {
//                    f[j] = (fL[j] + fR[j] - maxEigenValue * (WR[j] - WL[j])) / 2;
//                }

                double c2 = 1000;
                double[] Wstar = new double[nEqs];
                // pressure
                double VnL = 0;
                double VnR = 0;
                for (int d = 0; d < dim; ++d) {
                    VnL += WL[d + 1] * n[d];
                    VnR += WR[d + 1] * n[d];
                }
                double beta = (VnL + VnR) / 2;
                double alfa = Math.sqrt(beta * beta + 4 * c2);
                double dp = WR[0] - WL[0];
                double du = VnR - VnL;
                double us = (VnL + VnR) / 2 - dp / alfa - beta / (2 * alfa) * du;
                double ps = (WL[0] + WR[0]) / 2 - c2 / alfa * du + beta / (2 * alfa) * dp;
                
                // tangential velocity
                double[] vt = new double[dim];
                if (us > 0) {
                    for (int d = 0; d < dim; d++) {
                        vt[d] = WL[d + 1] - VnL * n[d];
                    }
                } else {
                    for (int d = 0; d < dim; d++) {
                        vt[d] = WR[d + 1] - VnR * n[d];
                    }
                }
                
                Wstar[0] = ps;
                for (int d = 0; d < dim; d++) {
                    Wstar[d + 1] = us * n[d] + vt[d];
                }
                f = convectiveFlux(Wstar, n, elem);
                break;
        }

        return f;
    }

    @Override
    public double[] convectiveFlux(double[] W, double[] n, ElementData elem) {
        double[] f = new double[nEqs];

        double V = .0;
        for (int d = 0; d < dim; ++d) {
            V += W[d + 1] * n[d];
        }

        f[0] = V;
        for (int d = 0; d < dim; ++d) {
            f[d + 1] = W[d + 1] * V + W[0] * n[d];
        }
        return f;
    }

    @Override
    public boolean isEquationsJacobian() {
        return false;
    }

    @Override
    public double[] diffusiveFlux(double[] W, double[] dW, double[] n, ElementData elem) {

        double[] velocityJac = new double[dim * dim];
        for (int d = 0; d < dim; ++d) {
            for (int f = 0; f < dim; ++f) {
                velocityJac[dim * d + f] = dW[f * nEqs + d + 1];
            }
        }

        // stress tensor calculation
        double[] stress = new double[dim * dim];
        for (int d = 0; d < dim; ++d) {
            for (int f = 0; f < dim; ++f) {
                stress[dim * d + f] = (velocityJac[dim * d + f] + velocityJac[dim * f + d]) / 2;
            }
        }

        double[] flux = new double[nEqs];
        flux[0] = 0;
        for (int d = 0; d < dim; ++d) {
            for (int f = 0; f < dim; ++f) {
                flux[f + 1] += stress[dim * d + f] * n[d] / Re;
            }
        }
        return flux;
    }

    @Override
    public double[] numericalDiffusiveFlux(double Wc[], double dWc[], double[] n, int TT, ElementData elem) {
        double[] flux = diffusiveFlux(Wc, dWc, n, elem);
        if (TT < -1) {
            Arrays.fill(flux, .0);
        }
        return flux;
    }

    @Override
    public boolean isSourcePresent() {
        return false;
    }

    @Override
    public double[] sourceTerm(double[] W, double[] dW, ElementData elem) { // zdrojovy clen
        throw new UnsupportedOperationException("source is not present");
    }

    @Override
    public double[] boundaryValue(double[] WL, double[] n, int TT, ElementData elem) {
        double[] WR = new double[nEqs];

        switch (TT) {
            case (BoundaryType.WALL):
                if (isDiffusive) {
                    double[] u = elem.meshVelocity;
                    WR[0] = WL[0];
                    for (int d = 0; d < dim; ++d) {
                        WR[d + 1] = WR[0] * u[d];
                    }
                } else {
                    WR = Arrays.copyOf(WL, nEqs);
                    double nu = 0;
                    for (int d = 0; d < dim; ++d) {
                        nu += WL[d + 1] * n[d];
                    }
                    for (int d = 0; d < dim; ++d) { //tangent to wall
                        WR[d + 1] = WL[d + 1] + n[d] * nu;
                    }
                }
                break;
            case (BoundaryType.INVISCID_WALL):
                WR = Arrays.copyOf(WL, nEqs);
                double nu = 0;
                for (int d = 0; d < dim; ++d) {
                    nu += WL[d + 1] * n[d];
                }
                for (int d = 0; d < dim; ++d) { //tangent to wall
                    WR[d + 1] = WL[d + 1] + n[d] * nu;
                }
                break;

            case (BoundaryType.INLET):
                WR[0] = WL[0];
                for (int d = 0; d < dim; ++d) {
                    WR[d + 1] = VIn[d];
                }
                break;

            case (BoundaryType.OUTLET):
                WR[0] = pOut;
                for (int d = 0; d < dim; ++d) {
                    WR[d + 1] = WL[d + 1];
                }
                break;
        }
        return WR;
    }

    @Override
    public double pressure(double[] W) {
        return 1;
    }

    @Override
    public double maxEigenvalue(double[] W, ElementData elem) {
        limite(W);
        double u = Math.sqrt(W[1] * W[1] + W[2] * W[2]);
        return u;
    }

    @Override
    public double[] convectiveFluxJacobian(double[] W, double[] n, ElementData elemData) {
        throw new UnsupportedOperationException("operation not supported");
    }

    @Override
    public double[] diffusiveFluxJacobian(double[] W, double[] dW, double n[], ElementData elemData) {
        throw new UnsupportedOperationException("operation not supported");
    }

    @Override
    public double[] sourceTermJacobian(double[] W, double[] dW, ElementData elemData) {
        throw new UnsupportedOperationException("operation not supported");
    }

    void limite(double[] W) {
    }

    @Override
    public boolean isIPFace(int TT) {
        return TT == BoundaryType.WALL;
    }

    @Override
    public void saveReferenceValues(String filePath) throws IOException {
        FlowProProperties output = new FlowProProperties();

        output.setProperty("l", Double.toString(lRef));
        output.setProperty("v", Double.toString(velocityRef));
        output.setProperty("t", Double.toString(tRef));

        output.store(new FileOutputStream(filePath), null);
    }

    @Override
    public double[] getReferenceValues() {
        return new double[]{lRef, velocityRef, tRef};
    }

    @Override
    public double[] getResults(double[] W, double[] dW, double[] X, String name) {
        switch (name) {
            case "pressure":
                return new double[]{W[0]};

            case "xVelocity":
                return new double[]{velocityRef * W[1]};

            case "yVelocity":
                if (dim > 1) {
                    return new double[]{velocityRef * W[2]};
                } else {
                    throw new UnsupportedOperationException("undefined value" + name);
                }

            case "zVelocity":
                if (dim > 3) {
                    return new double[]{velocityRef * W[3]};
                } else {
                    throw new UnsupportedOperationException("undefined value" + name);
                }

            case "velocity":
                double[] velocity = new double[dim];
                for (int i = 0; i < dim; i++) {
                    velocity[i] = velocityRef * W[i + 1];
                }
                return velocity;

            case "div":
                double[] div = new double[1];
                for (int i = 0; i < dim; i++) {
                    div[0] = dW[nEqs * i + i + 1];
                }
                return div;

            default:
                throw new UnsupportedOperationException("undefined value " + name);
        }
    }
}
