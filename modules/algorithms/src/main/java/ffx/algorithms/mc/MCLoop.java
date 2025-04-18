// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// ******************************************************************************
package ffx.algorithms.mc;

import static ffx.utilities.Constants.R;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.random;

import ffx.algorithms.dynamics.thermostats.Thermostat;
import ffx.algorithms.optimize.Minimize;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.utils.Loop;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

/**
 * MCLoop class.
 *
 * @author Armin Avdic
 */
public class MCLoop implements MonteCarloListener {

  private static final Logger logger = Logger.getLogger(MCLoop.class.getName());
  /** The MolecularAssembly. */
  private final MolecularAssembly molecularAssembly;
  /** The MD thermostat. */
  private final Thermostat thermostat;
  /**
   * First residue number of loop.
   */
  private final int firstResidue;
  /**
   * Last residue number of loop.
   */
  private final int endResidue;
  /** Everyone's favorite. */
  private final Random random = new Random();
  /** The ForceFieldEnergy object being used by MD. */
  private final ForceFieldEnergy forceFieldEnergy;
  /** KIC generation of loop solutions. See doi:10.1002/jcc.10416 */
  final Loop loop;
  /** The current MD step. */
  private int stepCount = 0;
  /** Number of simulation steps between MC move attempts. */
  private final int mcStepFrequency;
  /** Number of accepted MD moves. */
  private int numMovesAccepted;
  /** Number of KIC iterations per MC move. */
  private int iterations;
  /** List of active atoms. */
  private Atom[] atoms;

  private boolean skipAlgorithm = false;

  /**
   * Construct a Monte-Carlo loop switching mechanism.
   *
   * @param molecularAssembly the molecular assembly
   * @param mcStepFrequency number of MD steps between switch attempts
   * @param thermostat the MD thermostat
   * @param firstResidue first residue number of loop
   * @param endResidue last residue number of loop
   */
  MCLoop(MolecularAssembly molecularAssembly, int mcStepFrequency, Thermostat thermostat,
      int firstResidue, int endResidue) {
    numMovesAccepted = 0;

    this.molecularAssembly = molecularAssembly;
    this.atoms = molecularAssembly.getAtomArray();
    this.forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    this.mcStepFrequency = (mcStepFrequency == 0) ? Integer.MAX_VALUE : mcStepFrequency;
    this.thermostat = thermostat;
    /* Energy of the system at initialization. */
    double systemReferenceEnergy = molecularAssembly.getPotentialEnergy().energy(false, true);
    this.firstResidue = firstResidue;
    this.endResidue = endResidue;
    this.iterations = 1;

    if ((endResidue - firstResidue) < 3) {
      logger.info("MCLoop requires at least 3 residues. First and last residues are anchors.");
      skipAlgorithm = true;
    }
    String sb =
        " Running MCLoop:\n" + format("     mcStepFrequency: %4d\n", mcStepFrequency) + format(
            "     referenceEnergy: %7.2f\n", systemReferenceEnergy);
    logger.info(sb);

    loop = new Loop(molecularAssembly);
  }

  /**
   * Get the current MC acceptance rate.
   *
   * @return the acceptance rate.
   */
  public double getAcceptanceRate() {
    // Intentional integer division.
    int numTries = stepCount / mcStepFrequency;
    return (double) numMovesAccepted / numTries;
  }

  /**
   * {@inheritDoc}
   *
   * <p>The primary driver. Called by the MD engine at each dynamics step.
   */
  @Override
  public boolean mcUpdate(double temperature) {

    stepCount++;
    if (skipAlgorithm) {
      return false;
    }
    // Decide on the type of step to be taken.
    if ((stepCount % mcStepFrequency != 0)) {
      // Not yet time for an MC step, return to MD.
      return false;
    }

    atoms = molecularAssembly.getAtomArray();

    // Randomly choose a target sub portion of loop to KIC.
    int midResidue;
    midResidue = ThreadLocalRandom.current().nextInt(firstResidue + 1, endResidue);

    List<double[]> loopSolutions;
    loopSolutions = loop.generateLoops(midResidue - 1, midResidue + 1);

    for (int i = 1; i < iterations; i++) {
      // pick random subLoop
      midResidue = ThreadLocalRandom.current().nextInt(firstResidue + 1, endResidue);
      // pick random solution
      if (!loopSolutions.isEmpty()) {
        List<double[]> tempLoops = loop.generateLoops(midResidue - 1, midResidue + 1,
            loopSolutions.get(random.nextInt(loopSolutions.size())));
        loopSolutions.addAll(tempLoops);
      } else {
        loopSolutions = loop.generateLoops(midResidue - 1, midResidue + 1);
      }
    }
    int numLoopsFound = loopSolutions.size();
    // Check whether KIC found alternative loops
    if (numLoopsFound <= 1) {
      return false;
    }

    // Perform the MC move.
    return tryLoopStep(loopSolutions);
  }

  /**
   * Setter for the field <code>iterations</code>.
   *
   * @param iterations The number of KIC iterations per MC move.
   */
  public void setIterations(int iterations) {
    this.iterations = iterations;
  }

  /**
   * Perform a loop MC move.
   *
   * @param loopSolutions A list of loop solutions.
   * @return accept/reject
   */
  private boolean tryLoopStep(List<double[]> loopSolutions) {

    // Choose from the list of available loops and save current coordinates
    double[] newCoords = loopSolutions.get(random.nextInt(loopSolutions.size()));
    // Storage for coordinates before MC move.
    double[] oldCoords = storeActiveCoordinates();

    // Optimize the system.
    Minimize minimize1 = new Minimize(null, forceFieldEnergy, null);
    minimize1.minimize();
    double originalLoopEnergy = forceFieldEnergy.energy(false, true);

    // Perform move and analysis of chosen loop.
    performLoopMove(newCoords);
    Minimize minimize2 = new Minimize(null, forceFieldEnergy, null);
    minimize2.minimize();
    double newLoopEnergy = forceFieldEnergy.energy(false, true);

    double temperature = thermostat.getCurrentTemperature();
    double kT = R * temperature;
    // Test the MC criterion for a loop move.
    double dE = Math.abs(originalLoopEnergy - newLoopEnergy);
    if (newLoopEnergy < originalLoopEnergy) {
      dE = -dE;
    }

    StringBuilder sb = new StringBuilder();
    sb.append(" Assessing possible MC loop move step:\n");
    sb.append(format("     original loop: %16.8f\n", originalLoopEnergy));
    sb.append(format("     possible loop: %16.8f\n", newLoopEnergy));
    sb.append("     -----\n");

    // Test Monte-Carlo criterion.
    if (dE < 0) {
      sb.append("     Accepted!");
      logger.info(sb.toString());
      numMovesAccepted++;
      return true;
    }
    double criterion = exp(-dE / kT);

    double metropolis = random();
    sb.append(format("     criterion:  %9.4f\n", criterion));
    sb.append(format("     rng:        %9.4f\n", metropolis));
    if ((metropolis < criterion)) {
      sb.append("     Accepted!");
      logger.info(sb.toString());
      numMovesAccepted++;
      return true;
    }
    sb.append("     Denied.");
    logger.info(sb.toString());

    // Move was rejected, undo the loop move
    performLoopMove(oldCoords);
    return false;
  }

  /**
   * Perform the requested coordinate move
   *
   * @param newCoordinates THe new coordinates.
   */
  private void performLoopMove(double[] newCoordinates) {
    int index = 0;
    for (Atom a : atoms) {
      double x = newCoordinates[index++];
      double y = newCoordinates[index++];
      double z = newCoordinates[index++];
      a.moveTo(x, y, z);
    }
  }

  private double[] storeActiveCoordinates() {
    double[] coords = new double[atoms.length * 3];
    int index = 0;
    for (Atom a : atoms) {
      coords[index++] = a.getX();
      coords[index++] = a.getY();
      coords[index++] = a.getZ();
    }
    return coords;
  }

}
