package flowpro.user.equation;

import flowpro.api.ElementData;
import flowpro.api.Equation;
import flowpro.api.FlowProProperties;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

public class IncompressibleNavierStokesSpalartAllmaras implements Equation {

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
    protected double velocityRef;
    protected double tRef;
    protected double rho;
    protected double eta;

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
    
    double vtIn; // turbulence intensity at the inlet
    
    // parameters of turbulence model
    double sigma;
    double cb1;
    double cb2;
    double ka;
    double cw1;
    double cw2;
    double cw3;
    double cv1;
    double ct1;
    double ct2;
    double ct3;
    double ct4;
    double Prt;
    double C_prod;

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
        nEqs = dim + 2;

        gravityAcceleration = 9.81;

        isDiffusive = props.getBoolean("isFlowViscous");

        // reference values from inlet
        VIn = props.getDoubleArray("VIn");

        // density
        if (props.containsKey("density")){
            rho = props.getDouble("density");
        }
        
        if (props.containsKey("viscosity")){
            eta = props.getDouble("viscosity");
        }
        
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
        pRef = rho*velocityRef*velocityRef;
        
        for (int d = 0; d < dim; ++d) {
            VIn[d] /= velocityRef;
        }
        // outlet
        if (props.containsKey("pOut")) {
            pOut = props.getDouble("pOut");
        } else {
            throw new IOException("outlet boundary pressure condition is not specified");
        }

        if (props.containsKey("reynolds")){
            Re = props.getDouble("reynolds");
        } else {
            Re = rho*velocityRef*lRef/eta;
        }
        
        System.out.println("Reynolds number " + Re);
        
        // turbulence in inlet
        vtIn = props.getDouble("vtIn");
        
        // parameters of the turbulence model
        sigma = props.getDouble("sigma");
        cb1 = props.getDouble("cb1");
        cb2 = props.getDouble("cb2");
        ka = props.getDouble("ka");
        cw1 = cb1 / (ka * ka) + (1 + cb2) / sigma;
        cw2 = props.getDouble("cw2");
        cw3 = props.getDouble("cw3");
        cv1 = props.getDouble("cv1");
        ct1 = props.getDouble("ct1");
        ct2 = props.getDouble("ct2");
        ct3 = props.getDouble("ct3");
        ct4 = props.getDouble("ct4");
        Prt = props.getDouble("Prt");
        C_prod = props.getDouble("C_prod");
    }

    @Override
    public void setState(double dt, double t) {
    }

    @Override
    public double[] constInitCondition() {
        if (dim == 2) {
            return new double[]{pOut, VIn[0], VIn[1], vtIn};
        } else {
            return new double[]{pOut, VIn[0], VIn[1], VIn[2], vtIn};
        }
    }

    //  nevazky tok stenou _____________________________________________________
    @Override
    public double[] numericalConvectiveFlux(double WL[], double WR[], double[] n, int TT, ElementData elem) {
        double[] f = new double[nEqs];

        switch (TT) {
            case (BoundaryType.WALL):
            case (BoundaryType.INVISCID_WALL):
                
                f[0] = 0;
                for (int d = 0; d < dim; d++) {
                    f[d + 1] = WR[0] * n[d];
                }
                
                // for ALE
                double V = 0;
                for (int d = 0; d < dim; d++) {
                    V += elem.meshVelocity[d] * n[d];
                }
                f[0] += V;
                for (int j = 1; j < nEqs; j++) {
                    f[j] += V * WR[j];
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
                double vts = (WL[dim + 1] + WR[dim + 1])/2;

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
                Wstar[dim + 1] = vts;
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
        f[dim + 1] = W[dim + 1] * V;
        
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

        double vt = max(0, W[dim + 1]);
        double[] vtDer = new double[dim];
        for (int d = 0; d < dim; ++d) {
            vtDer[d] = dW[d * nEqs + dim + 1];
        }

        double xi = vt;
        double fv1 = xi * xi * xi / (xi * xi * xi + cv1 * cv1 * cv1);
        double mut = rho * vt * fv1; // nejsem si jisty tim rho !!!!!!!!!!!
        
        double[] flux = new double[nEqs];
        flux[0] = 0;
        for (int d = 0; d < dim; d++) {
            for (int f = 0; f < dim; f++) {
                flux[f + 1] += (1 + mut) * stress[dim * d + f] * n[d] / Re;
            }
            flux[dim + 1] += 1 / (Re * sigma) * (1 + vt) * vtDer[d] * n[d];
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
        return true;
    }

    @Override
    public double[] sourceTerm(double[] W, double[] dW, ElementData elem) { // zdrojovy clen
        
        double[] velocityJac = new double[dim * dim];
        for (int d = 0; d < dim; ++d) {
            for (int f = 0; f < dim; ++f) {
                velocityJac[dim * d + f] = dW[f * nEqs + d + 1];
            }
        }

        double vt = max(0, W[dim + 1]);
        double vtDerMag = 0;
        for (int d = 0; d < dim; ++d) {
            double vtDer = dW[d * nEqs + dim + 1];
            vtDerMag += vtDer * vtDer;
        }

        // turbulence limit
        if (vt < 0) {
            vt = 0;
        }

        double D = elem.currentWallDistance;

        double xi = vt; // vt/v 
        double ft2 = 0; //ct3*Math.exp(-ct4*xi*xi);
        double ft1 = 0;
        double fv1 = xi * xi * xi / (xi * xi * xi + cv1 * cv1 * cv1);
        double fv2 = 1 - xi / (1 + xi * fv1);
        double Om = rotationMagnitude(velocityJac);
        double S = Om + 1 / Re * vt / (ka * ka * D * D) * fv2;
        double rt = vt / (Re * S * ka * ka * D * D);
        if (rt > 10) {
            rt = 10;
        }
        double g = rt + cw2 * (Math.pow(rt, 6.0) - rt);
        double fw = g * Math.pow((1 + Math.pow(cw3, 6.0)) / (Math.pow(g, 6.0) + Math.pow(cw3, 6.0)), 1.0 / 6);

        double[] source = new double[nEqs];
        source[dim + 1] = limitDestruction(1 / Re * cb2 * vtDerMag + cb1 * (1 - ft2) * S * vt - 1 / Re * (cw1 * fw - cb1 / (ka * ka) * ft2) * (vt / D) * (vt / D)); // + r*Re*ft1*dU*dU
        return source;
    }

    @Override
    public double[] boundaryValue(double[] WL, double[] n, int TT, ElementData elem) {
        double[] WR = new double[nEqs];

        switch (TT) {
            case (BoundaryType.WALL):
                if (isDiffusive) {
                    double[] u = elem.meshVelocity;
                    WR[0] = WL[0];
                    for (int d = 0; d < dim; d++) {
                        WR[d + 1] = u[d];
                    }
                    WR[dim + 1] = 0;
                } else {
                    WR[0] = WL[0];
                    double nu = 0;
                    for (int d = 0; d < dim; d++) {
                        nu += WL[d + 1] * n[d];
                    }
                    for (int d = 0; d < dim; ++d) { //tangent to wall
                        WR[d + 1] = WL[d + 1] - n[d] * nu;
                    }
                    WR[dim + 1] = WL[dim + 1];
                }
                break;
            case (BoundaryType.INVISCID_WALL):
                WR[0] = WL[0];
                double nu = 0;
                for (int d = 0; d < dim; d++) {
                    nu += WL[d + 1] * n[d];
                }
                for (int d = 0; d < dim; d++) { //tangent to wall
                    WR[d + 1] = WL[d + 1] - n[d] * nu;
                }
                WR[dim + 1] = WL[dim + 1];
                break;

            case (BoundaryType.INLET):
                WR[0] = WL[0];
                for (int d = 0; d < dim; ++d) {
                    WR[d + 1] = VIn[d];
                }
                WR[dim + 1] = vtIn;
                break;

            case (BoundaryType.OUTLET):
                WR[0] = pOut;
                for (int d = 0; d < dim; ++d) {
                    WR[d + 1] = WL[d + 1];
                }
                WR[dim + 1] = WL[dim + 1];
                break;
        }
        return WR;
    }

    @Override
    public double pressure(double[] W) {
        return W[0];
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
    public double[] combineShockSensors(double[] shock) {
        for (int m = 1; m < nEqs; m++) {
            shock[m] = shock[0]; // all shock sensors are acording divergence
        }
        return shock;
    }

    @Override
    public void saveReferenceValues(String filePath) throws IOException {
        FlowProProperties output = new FlowProProperties();

        output.setProperty("l", Double.toString(lRef));
        output.setProperty("v", Double.toString(velocityRef));
        output.setProperty("t", Double.toString(tRef));
        output.setProperty("rho", Double.toString(rho));
        output.setProperty("p", Double.toString(pRef));

        output.store(new FileOutputStream(filePath), null);
    }

    @Override
    public double[] getReferenceValues() {
        return new double[]{lRef, pRef, rho, velocityRef, tRef};
    }

    @Override
    public double[] getResults(double[] W, double[] dW, double[] X, String name) {
        switch (name) {
            case "pressure":
                return new double[]{pRef*W[0]};

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

            case "mut":
                double xi = max(W[dim + 1],0); // vt/v
                double fv1 = xi * xi * xi / (xi * xi * xi + cv1 * cv1 * cv1);
                return new double[]{W[dim + 1]*fv1};    
                
            default:
                throw new UnsupportedOperationException("undefined value " + name);
        }
    }
    
    public double max(double a, double b) {
        if (a > b) {
            return a;
        } else {
            return b;
        }
    }

    public double matrixMagnitude(double[] A) {
        double mag = 0;
        for (int i = 0; i < A.length; i++) {
            mag += A[i] * A[i];
        }
        return Math.sqrt(mag);
    }

    public double rotationMagnitude(double[] U) {
        double rotMag = 0;
        if (U.length == 4) {
            rotMag = Math.abs(U[1] - U[2]);
        } else if (U.length == 9) {
            rotMag = Math.sqrt((U[0*dim + 1] - U[1*dim + 0]) * (U[0*dim + 1] - U[1*dim + 0]) + (U[0*dim + 2] - U[2*dim + 0]) * (U[0*dim + 2] - U[2*dim + 0]) + (U[2*dim + 1] - U[1*dim + 2]) * (U[2*dim + 1] - U[1*dim + 2]));
        }
        return rotMag;
    }

    public double limitDestruction(double d) {
        if (d < -100) {
            d = -100;
        }
        return d;
    }
}
