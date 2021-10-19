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
	
	protected enum InletType {
		VELOCITY, PARABOLA, PARABOLA2
	}
	
	protected InletType inletType;
		
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
    // inlet boundary condition 
    protected final double pIn0 = 1; // static pressure
    protected final double rhoIn0 = 1; // density
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
		
		if (props.containsKey("inletType")) {
			String type = props.getString("inletType").toUpperCase();
			
            if (InletType.VELOCITY.name().toLowerCase().equals(type)) {
				inletType = InletType.VELOCITY;
			} else if (InletType.PARABOLA.name().toLowerCase().equals(type)) {
				inletType = InletType.PARABOLA;
			} else if (InletType.PARABOLA2.name().toLowerCase().equals(type)) {
				inletType = InletType.PARABOLA2;
			} else {
				throw new IOException("unknown inlet type \'" + type + "\'");
			}
        } else {
            inletType = InletType.VELOCITY;
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
                //System.out.println(elem.meshVelocity[1]);
                double V = 0;
                for (int d = 0; d < dim; d++) {
                    V += elem.meshVelocity[d] * n[d];
                }
                f[0] += V;
                for (int d = 0; d < dim; d++) {
                    f[d + 1] += V * WR[d + 1];
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
////					f[j] = (fL[j] + fR[j]) / 2;
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
        double[] stress = viscousStressTensor(W, dW);
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
                    //System.out.println(u[0] + " " + u[1]);
                    WR[0] = WL[0];
					System.arraycopy(u, 0, WR, 1, dim);
                } else {
                    WR[0] = WL[0];
                    double nu = 0;
                    for (int d = 0; d < dim; d++) {
                        nu += WL[d + 1] * n[d];
                    }
                    for (int d = 0; d < dim; ++d) { //tangent to wall
                        WR[d + 1] = WL[d + 1] - n[d] * nu;
                    }
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
                break;

            case (BoundaryType.INLET):
                WR[0] = WL[0];				
                for (int d = 0; d < dim; ++d) {
                    WR[d + 1] = VIn[d];
                }
//				double height = 0.41;
//				double y = elem.currentX[1];
//				WR[1] = 1.5 * y * (height-y) / (height*height/4);
//				if (inletType == InletType.VELOCITY) {
//					for (int d = 0; d < dim; ++d) {
//						WR[d + 1] = VIn[d];
//					}
//				} else {
//					double height = 0.41;
//					double y = elem.currentX[1];
//					WR[1] = 1.5 * y * (height-y) / (height*height/4);
//				}
//				
//				double t = elem.currentT * tRef;
//				if (inletType == InletType.PARABOLA2 && t < 2.0) {
//					WR[1] *= (1 - Math.cos(Math.PI / 2 * t)) / 2;
//				}
				
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
        return W[0];
    }
	
	protected double[] viscousStressTensor(double[] W, double[] dW) {
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
                stress[dim * d + f] = (velocityJac[dim * d + f] + velocityJac[dim * f + d]);
            }
        }
		
		return stress;
	}

    @Override
	public double[] stressVector(double[] W, double[] dW, double[] normal) {	
		double p = pressure(W);
		double[] viscousStress = viscousStressTensor(W, dW);
		double[] stressVector = new double[dim];
		for (int d = 0; d < dim; ++d) {
			stressVector[d] -= p * normal[d];
			for (int f = 0; f < dim; ++f) {
				stressVector[d] += 1 / Re * viscousStress[d * dim + f] * normal[f];
			}
		}		
		
		return stressVector;
	}
    
    @Override
    public double maxEigenvalue(double[] W, ElementData elem) {
//        limite(W);
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

//    void limite(double[] W) {
//    }

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
                    div[0] += dW[nEqs * i + i + 1];
                }
                return div;
            
            case "vorticity":
                if (dim == 2) {
                    double dvdx = dW[2];
                    double dudy = dW[nEqs + 1];
                    return new double[] {velocityRef / lRef * (dvdx - dudy)};                
                } else {
                    throw new UnsupportedOperationException("quantity \"" + name
                            + "\" is only available in two dimensions");
                }
            default:
                throw new UnsupportedOperationException("undefined value " + name);
        }
    }
}
