package population_model;

import java.awt.Color;
import java.util.Vector;

import javax.vecmath.Vector3d;

import processing.core.PGraphics3D;

import bsim.BSim;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimP3DDrawer;
import bsim.particle.BSimBacterium;
import bsim.ode.BSimOdeSolver;
import bsim.ode.BSimOdeSystem;
import bsim.export.BSimLogger;
import bsim.export.BSimMovExporter;
import bsim.export.BSimPngExporter;


public class BSimEColi {
  public static void main(String[] args) {

    BSim sim = new BSim();
    sim.setDt(0.01);
    sim.setSimulationTime(10.0);
    sim.setTimeFormat("0.00");
    sim.setBound(100, 100, 100);
    sim.setSolid(true, true, true);

    // Rate Constants
    double kcat_Nap = 2, kd_Nap = 0.6, km_Nap = 1;
    double kcat_Nrf = 2, kd_Nrf = 0.6, km_Nrf = 1;

    double k1 = kcat_Nap / (km_Nap - kd_Nap);
    double k_1 = kd_Nap * k1;
    double k3 = kcat_Nrf / (km_Nrf - kd_Nrf);
    double k_3 = kd_Nrf * k3;
    double k2 = kcat_Nap, k4 = kcat_Nrf;

    class EColiODEBacterium extends BSimBacterium {
      protected MyODE ode;
      protected double[] x, xNew;
      public EColiODEBacterium(BSim sim, Vector3d position) {
        super(sim, position);
        ode = new MyODE();
        x = ode.getICs();
      }

      @Override
      public void action() {
        super.action();
        xNew = BSimOdeSolver.rungeKutta45(ode, sim.getTime(), x, sim.getDt());
        x = xNew;
      }

      class MyODE implements BSimOdeSystem {
        private int numEq = 7;
        public double[] derivativeSystem(double t, double[] x) {
          double[] dx = new double[numEq];
          // Set of ODES describing enzyme kinetics
          dx[0] = k_1 * x[2] - k1 * x[1] * x[0];
          dx[1] = k_1 * x[2] - k1 * x[1] * x[0] + k2 * x[2];
          dx[2] = k1 * x[1] * x[0] - k_1 * x[2] - k2 * x[2];
          dx[3] = k2 * x[2] + k_3 * x[5] - k3 * x[4] * x[3];
          dx[4] = k_3 * x[5] - k3 * x[4] * x[3] + k4 * x[5];
          dx[5] = k3 * x[3] * x[4] - k_3 * x[5] - k4 * x[5];
          dx[6] = k4 * x[5];

          return dx;
        }
        public double[] getICs() {
          double[] ics = new double[numEq];
          // Set up initial conditions
          ics[0] = 10; //10e-2;  // NO3
          ics[1] = 5; //3e-2;   // Nap
          ics[2] = 0.0;    // NapNO3
          ics[3] = 10;//5e-3;   // NO2
          ics[4] = 5;//3e-2;   // Nrf
          ics[5] = 0.0;    // NrfNO2
          ics[6] = 0.0;    // NH4
          return ics;
        }
        public int getNumEq() { return numEq; }
      }
      public double[] getODEValue() { return x; }
    }

    final Vector<EColiODEBacterium> bacteria = new Vector<EColiODEBacterium>();


    while(bacteria.size() < 1000) {
			EColiODEBacterium b = new EColiODEBacterium(sim,
			                      new Vector3d(Math.random()*sim.getBound().x,
						                             Math.random()*sim.getBound().y,
						                             Math.random()*sim.getBound().z));
			if(!b.intersection(bacteria)) bacteria.add(b);
		}

    sim.setTicker(new BSimTicker() {
      @Override
      public void tick() {
        for(EColiODEBacterium b : bacteria) {
          b.action();
          b.updatePosition();
        }
      }
    });

    BSimP3DDrawer drawer = new BSimP3DDrawer(sim, 1000,1000) {
			@Override
			public void scene(PGraphics3D p3d) {
				for(EColiODEBacterium b : bacteria) {
          int G = (int) (10 * b.getODEValue()[6]);
          int R = 150;
          int B = 150;
          draw(b, new Color(R, G, B));
				}
			}
		};
		sim.setDrawer(drawer);

    String resultsDir = BSimUtils.generateDirectoryPath("./results/");

    BSimLogger logger = new BSimLogger(sim, resultsDir + "logger.csv") {
      @Override
      public void before() {
        super.before();
        String buffer = new String();
        for(int i = 0; i < 7; i++) {
          buffer = buffer + ",Equation_" + i + "_stat";
        }
        write("Time Seconds" + buffer);
      }

			@Override
			public void during() {
        String o = sim.getFormattedTime();
        String buffer = new String();
        //for (EColiODEBacterium b : bacteria) {
          for (int i = 0; i < 7; i++) {
            buffer = buffer + "," + bacteria.get(0).getODEValue()[i]; // Writes out value of NO3
          }
        //}
        write(o + buffer);
			}
		};
		logger.setDt(0.01);
		sim.addExporter(logger);

    BSimMovExporter movExporter = new BSimMovExporter(sim, drawer, resultsDir + "BSim.mov");
		movExporter.setDt(0.03);
		sim.addExporter(movExporter);

    BSimPngExporter pngExporter = new BSimPngExporter(sim, drawer, resultsDir);
		pngExporter.setDt(0.5);
		sim.addExporter(pngExporter);

    //sim.preview();
    sim.export();
  }
}
