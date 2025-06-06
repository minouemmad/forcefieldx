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
package ffx.algorithms.dynamics;

import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;
import static org.apache.commons.math3.util.FastMath.exp;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import java.io.IOException;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The ReplicaExchange implements temperature and lambda replica exchange methods.
 *
 * @author Timothy D. Fenn and Michael J. Schnieders
 * @since 1.0
 */
public class ReplicaExchange implements Terminatable {

  private static final Logger logger = Logger.getLogger(ReplicaExchange.class.getName());
  private final int nReplicas;
  private final Random random;
  /** Parallel Java world communicator. */
  private final Comm world;
  /** Rank of this process. */
  private final int rank;
  /**
   * The parameters array stores communicated parameters for each process (i.e. each RepEx system).
   * Currently, the array is of size [number of Processes][2].
   */
  private final double[][] parameters;
  /**
   * Each parameter array is wrapped inside a Parallel Java DoubleBuf for the All-Gather
   * communication calls.
   */
  private final DoubleBuf[] parametersBuf;
  private final MolecularDynamics replica;
  private boolean done = false;
  private boolean terminate = false;
  private final double[] myParameters;
  private final DoubleBuf myParametersBuf;
  private final int[] temp2Rank;
  private final int[] rank2Temp;
  private double[] temperatures;
  private final int[] tempAcceptedCount;
  private final int[] rankAcceptedCount;
  private final int[] tempTrialCount;
  boolean monteCarlo;

  /**
   * ReplicaExchange constructor.
   *
   * @param molecularDynamics The MolecularDynamics instance.
   * @param listener A listener for algorithm events.
   * @param temperature The temperature (K).
   * @param exponent a double to set temperature ladder.
   */
  public ReplicaExchange(MolecularDynamics molecularDynamics, AlgorithmListener listener,
      double temperature, double exponent, boolean monteCarlo) {

    this.replica = molecularDynamics;
    this.monteCarlo = monteCarlo;

    // MolecularAssembly[] molecularAssemblies = molecularDynamics.getAssemblies();
    // CompositeConfiguration properties = molecularAssemblies[0].getProperties()

    // Set up the Replica Exchange communication variables for Parallel Java communication between
    // nodes.
    world = Comm.world();

    // Number of processes is equal to the number of replicas.
    int numProc = world.size();
    rank = world.rank();

    nReplicas = numProc;
    temperatures = new double[nReplicas];
    temp2Rank = new int[nReplicas];
    rank2Temp = new int[nReplicas];
    tempAcceptedCount = new int[nReplicas];
    rankAcceptedCount = new int[nReplicas];
    tempTrialCount = new int[nReplicas];

    setExponentialTemperatureLadder(temperature, exponent);

    random = new Random();
    random.setSeed(0);

    // Create arrays to store the parameters of all processes.
    parameters = new double[nReplicas][2];
    parametersBuf = new DoubleBuf[nReplicas];
    for (int i = 0; i < nReplicas; i++) {
      parametersBuf[i] = DoubleBuf.buffer(parameters[i]);
    }

    // A convenience reference to the parameters of this process are updated
    // during communication calls.
    myParameters = parameters[rank];
    myParametersBuf = parametersBuf[rank];
  }

  /**
   * Sample.
   *
   * @param cycles The number of cycles to sample.
   * @param nSteps The number of steps per cycle.
   * @param timeStep The time step (fsec).
   * @param printInterval The interval (in steps) to print current status.
   * @param saveInterval The interval (in steps) to save a snapshot.
   */
  public void sample(int cycles, long nSteps, double timeStep, double printInterval,
      double saveInterval) {
    done = false;
    terminate = false;
    for (int i = 0; i < cycles; i++) {
      // Check for termination request.
      if (terminate) {
        done = true;
        break;
      }
      dynamic(nSteps, timeStep, printInterval, saveInterval);
      logger.info(String.format(" Applying exchange condition for cycle %d.", i));
      exchange(i);
    }
  }

  /**
   * setExponentialTemperatureLadder.
   *
   * @param lowTemperature a double.
   * @param exponent a double.
   */
  public void setExponentialTemperatureLadder(double lowTemperature, double exponent) {
    for (int i = 0; i < nReplicas; i++) {
      temperatures[i] = lowTemperature * exp(exponent * i);
      temp2Rank[i] = i;
      rank2Temp[i] = i;
    }
  }

  /**
   * Setter for the field <code>temperatures</code>.
   *
   * @param temperatures an array of {@link double} objects.
   */
  public void setTemperatures(double[] temperatures) {
    assert (temperatures.length == nReplicas);
    this.temperatures = temperatures;
  }

  /**
   * {@inheritDoc}
   *
   * <p>This should be implemented as a blocking interrupt; when the method returns the <code>
   * Terminatable</code> algorithm has reached a clean termination point. For example, between
   * minimize or molecular dynamics steps.
   */
  @Override
  public void terminate() {
    terminate = true;
    while (!done) {
      synchronized (this) {
        try {
          wait(1);
        } catch (InterruptedException e) {
          logger.log(Level.WARNING, "Exception terminating replica exchange.\n", e);
        }
      }
    }
  }

  /**
   * All processes complete the exchanges identically given the same Random number seed.
   *
   * @param cycle The current cycle.
   */
  private void exchange(int cycle) {
    int start;
    int increment;
    if (monteCarlo) {
      // 1 M.C. trial per temperature (except for those at the ends of the ladder).
      start = cycle % 2;
      increment = 2;
    } else {
      // 2 M.C. trials per temperature (except for those at the ends of the ladder).
      start = 0;
      increment = 1;
    }
    // Loop over temperatures
    for (int temperature = start; temperature < nReplicas - 1; temperature += increment) {

      // Ranks for temperatures A and B
      int rankA = temp2Rank[temperature];
      int rankB = temp2Rank[temperature + 1];

      // Load temperature, beta and energy for each rank.
      double tempA = parameters[rankA][0];
      double tempB = parameters[rankB][0];
      double betaA = KCAL_TO_GRAM_ANG2_PER_PS2 / (tempA * kB);
      double betaB = KCAL_TO_GRAM_ANG2_PER_PS2 / (tempB * kB);
      double energyA = parameters[rankA][1];
      double energyB = parameters[rankB][1];

      // Compute the change in energy over kT (E/kT) for the Metropolis criteria.
      double deltaE = (energyA - energyB) * (betaB - betaA);
      // If the Metropolis criteria is satisfied, do the switch.

      //Count the number of trials for each temp
      tempTrialCount[temperature]++;
      if (deltaE < 0.0 || random.nextDouble() < exp(-deltaE)) {
        tempAcceptedCount[temperature]++;
        double tempAcceptance =
            tempAcceptedCount[temperature] * 100.0 / (tempTrialCount[temperature]);

        double rankAcceptance;
        if (tempA < tempB) {
          rankAcceptedCount[rankA]++;
          rankAcceptance = rankAcceptedCount[rankA] * 100.0 / (tempTrialCount[temperature]);
        } else {
          rankAcceptedCount[rankB]++;
          rankAcceptance = rankAcceptedCount[rankB] * 100.0 / (tempTrialCount[temperature + 1]);
        }

        // Swap temperature and energy values.
        parameters[rankA][0] = tempB;
        parameters[rankB][0] = tempA;

        parameters[rankA][1] = energyB;
        parameters[rankB][1] = energyA;

        // Map temperatures to process ranks.
        temp2Rank[temperature] = rankB;
        temp2Rank[temperature + 1] = rankA;

        // Map ranks to temperatures.
        rank2Temp[rankA] = temperature + 1;
        rank2Temp[rankB] = temperature;

        logger.info(String.format(
            " RepEx accepted (%5.1f%%) (%5.1f%%) for %6.2f (%d) and %6.2f (%d) for dE=%10.4f.",
            tempAcceptance, rankAcceptance, tempA, rankA, tempB, rankB, deltaE));
      } else {
        double tempAcceptance =
            tempAcceptedCount[temperature] * 100.0 / (tempTrialCount[temperature]);
        double rankAcceptance =
            rankAcceptedCount[temperature] * 100.0 / (tempTrialCount[temperature]);
        logger.info(String.format(
            " RepEx rejected (%5.1f%%) (f%5.1f%%) for %6.2f (%d) and %6.2f (%d) for dE=%10.4f.",
            tempAcceptance, rankAcceptance, tempA, rankA, tempB, rankB, deltaE));
      }
    }
  }

  /**
   * Blocking dynamic steps: when this method returns each replica has completed the requested number
   * of steps.
   *
   * @param nSteps the number of time steps.
   * @param timeStep the time step (fsec).
   * @param printInterval the number of steps between logging updates.
   * @param saveInterval the number of steps between saving snapshots.
   */
  private void dynamic(final long nSteps, final double timeStep, final double printInterval,
      final double saveInterval) {

    int i = rank2Temp[rank];

    // Start this processes MolecularDynamics instance sampling.
    boolean initVelocities = true;
    replica.dynamic(nSteps, timeStep, printInterval, saveInterval, temperatures[i], initVelocities,
        null);

    // Update this ranks' parameter array to be consistent with the dynamics.
    myParameters[0] = temperatures[i];
    myParameters[1] = replica.state.getPotentialEnergy();

    // Gather all parameters from the other processes.
    try {
      world.allGather(myParametersBuf, parametersBuf);
    } catch (IOException ex) {
      String message = " Replica Exchange allGather failed.";
      logger.log(Level.SEVERE, message, ex);
    }
  }
}
