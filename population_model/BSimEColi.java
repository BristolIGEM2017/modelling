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
import bsim.export.BSimMovExporter;
import bsim.export.BSimPngExporter;


public class BSimEColi {
  public static void main(String[] args) {

    BSim sim = new BSim();
    sim.setDt(0.01);
    sim.setSimulationTime(10.0);
    sim.setTimeFormat("0.00");
    sim.setBound(100, 100, 100);

    class EColiODEBacterium extends BSimBacterium {
      protected MyODE ode;
      protected double[] y, yNew;
      public EColiODEBacterium(BSim sim, Vector3d position) {
        super(sim, position);
        ode = new MyODE();
        y = ode.getICs();
      }

      @Override
      public void action() {
        super.action();
        yNew = BSimOdeSolver.rungeKutta45(ode, sim.getTime(), y, sim.getDt());
        y = yNew;
      }

      class MyODE implements BSimOdeSystem {
        private int numEq = 1;
        public double[] derivativeSystem(double x, double[] y) {
          double[] dy = new double[numEq];
          dy[0] = -Math.sin(x) + Math.cos(x);
          return dy;
        }
        public double[] getICs() {
          double[] ics = new double[numEq];
          ics[0] = 0.01;
          return ics;
        }
        public int getNumEq() {
          return numEq;
        }
      }
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

    BSimP3DDrawer drawer = new BSimP3DDrawer(sim, 800,600) {
			@Override
			public void scene(PGraphics3D p3d) {
				for(EColiODEBacterium b : bacteria) {
          int G = 150 + (int)(100.0 * (b.y[0] * 1.0));
          int R = 100 + (int)(100.0 * (b.y[0] * 1.0));
          int B = 50;
          if(R < 0) R = 0;
          if(G < 0) G = 0;
          draw(b, new Color(R, G, B));
				}
			}
		};
		sim.setDrawer(drawer);

    String resultsDir = BSimUtils.generateDirectoryPath("./results/");

    BSimMovExporter movExporter = new BSimMovExporter(sim, drawer, resultsDir + "BSim.mov");
		movExporter.setDt(0.03);
		sim.addExporter(movExporter);

    BSimPngExporter pngExporter = new BSimPngExporter(sim, drawer, resultsDir);
		pngExporter.setDt(0.5);
		sim.addExporter(pngExporter);

    sim.preview();
    //sim.export();
  }
}
