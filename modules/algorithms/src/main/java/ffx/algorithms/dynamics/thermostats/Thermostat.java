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
package ffx.algorithms.dynamics.thermostats;

import ffx.numerics.Constraint;
import ffx.numerics.Potential.VARIABLE_TYPE;
import ffx.potential.SystemState;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The abstract Thermostat class implements methods common to all thermostats for initializing
 * velocities from a Maxwell-Boltzmann distribution and computing the instantaneous temperature.
 * Abstract methods are declared for half-step and full-step modification of velocities for
 * thermostat implementations.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class Thermostat {

  private static final Logger logger = Logger.getLogger(Thermostat.class.getName());
  /**
   * The center of mass coordinates.
   */
  private final double[] centerOfMass = new double[3];
  /**
   * The linear momentum.
   */
  private final double[] linearMomentum = new double[3];
  /**
   * The angular momentum.
   */
  private final double[] angularMomentum = new double[3];
  /**
   * Number of degrees of freedom removed by constraints.
   */
  private final int constrainedDoF;
  /**
   * The identity of this Thermostat.
   */
  protected ThermostatEnum name;
  /**
   * The value of kT in kcal/mol at the target temperature.
   */
  protected double kT;
  /**
   * Number of degrees of freedom, which can be less than the number of variables. For example,
   * removing translational motion removes 3 degrees of freedom.
   */
  protected int degreesOfFreedom;
  /**
   * The molecular dynamics state to be used.
   */
  protected final SystemState state;
  /**
   * The type of each variable.
   */
  protected VARIABLE_TYPE[] type;
  /**
   * The random number generator that the Thermostat will use to initialize velocities.
   */
  protected Random random;
  /**
   * Any geometric constraints to apply during integration.
   */
  protected List<Constraint> constraints;
  /**
   * The target temperature that this thermostat should maintain.
   */
  double targetTemperature;
  /**
   * Flag to indicate that center of mass motion should be removed.
   */
  private boolean removeCenterOfMassMotion;
  /**
   * Reduce logging.
   */
  private boolean quiet = false;

  /**
   * Constructor for Thermostat.
   *
   * @param state             The molecular dynamics state to be used.
   * @param type              the VARIABLE_TYPE of each variable.
   * @param targetTemperature a double.
   */
  public Thermostat(SystemState state, VARIABLE_TYPE[] type, double targetTemperature) {
    this(state, type, targetTemperature, new ArrayList<>());
  }

  public Thermostat(SystemState state, VARIABLE_TYPE[] type, double targetTemperature, List<Constraint> constraints) {
    this.state = state;
    this.type = type;
    int n = state.getNumberOfVariables();
    assert (n > 3);
    assert (type.length == n);
    random = new Random();
    setTargetTemperature(targetTemperature);

    this.constraints = new ArrayList<>(constraints);
    // Not every type of constraint constrains just one DoF.
    // SETTLE constraints, for example, constrain three.
    double nConstrained = constraints.stream().mapToInt(Constraint::getNumDegreesFrozen).sum();

    double[] mass = state.getMass();
    int massConstrained = 0;
    for (int i = 0; i < n; i++) {
      if (mass[i] <= 0.0) {
        massConstrained++;
      }
    }

    if (massConstrained > 0 && nConstrained > 0) {
      logger.severe("Mass-constraints with other constraints are not supported.");
    }
    constrainedDoF = (int) max(nConstrained, massConstrained);

    removeCenterOfMassMotion = true;

    // Remove center of mass motion degrees of freedom.
    degreesOfFreedom = n - 3 - constrainedDoF;

    // Update the kinetic energy.
    computeKineticEnergy();
  }

  /**
   * Parse a string into a Thermostat enumeration.
   *
   * @param str Thermostat String.
   * @return An instance of the ThermostatEnum.
   */
  public static ThermostatEnum parseThermostat(String str) {
    try {
      return ThermostatEnum.valueOf(str.toUpperCase());
    } catch (Exception e) {
      logger.info(format(" Could not parse %s as a thermostat; defaulting to Berendsen.", str));
      return ThermostatEnum.BERENDSEN;
    }
  }

  /**
   * Compute the center of mass, linear momentum and angular momentum.
   *
   * @param remove If true, the center of mass motion will be removed.
   * @param print  If true, the center of mass and momenta will be printed.
   */
  public void centerOfMassMotion(boolean remove, boolean print) {
    double totalMass = 0.0;
    for (int i = 0; i < 3; i++) {
      centerOfMass[i] = 0.0;
      linearMomentum[i] = 0.0;
      angularMomentum[i] = 0.0;
    }

    int nVariables = state.getNumberOfVariables();
    double[] x = state.x();
    double[] v = state.v();
    double[] mass = state.getMass();

    int index = 0;
    while (index < nVariables) {
      double m = mass[index];
      if (m <= 0.0 || type[index] == VARIABLE_TYPE.OTHER) {
        index++;
        continue;
      }
      assert (type[index] == VARIABLE_TYPE.X);
      double xx = x[index];
      double vx = v[index++];
      assert (type[index] == VARIABLE_TYPE.Y);
      double yy = x[index];
      double vy = v[index++];
      assert (type[index] == VARIABLE_TYPE.Z);
      double zz = x[index];
      double vz = v[index++];
      totalMass += m;
      centerOfMass[0] += xx * m;
      centerOfMass[1] += yy * m;
      centerOfMass[2] += zz * m;
      linearMomentum[0] += vx * m;
      linearMomentum[1] += vy * m;
      linearMomentum[2] += vz * m;
      angularMomentum[0] += (yy * vz - zz * vy) * m;
      angularMomentum[1] += (zz * vx - xx * vz) * m;
      angularMomentum[2] += (xx * vy - yy * vx) * m;
    }

    angularMomentum[0] -= (centerOfMass[1] * linearMomentum[2] - centerOfMass[2] * linearMomentum[1]) / totalMass;
    angularMomentum[1] -= (centerOfMass[2] * linearMomentum[0] - centerOfMass[0] * linearMomentum[2]) / totalMass;
    angularMomentum[2] -= (centerOfMass[0] * linearMomentum[1] - centerOfMass[1] * linearMomentum[0]) / totalMass;
    centerOfMass[0] /= totalMass;
    centerOfMass[1] /= totalMass;
    centerOfMass[2] /= totalMass;
    linearMomentum[0] /= totalMass;
    linearMomentum[1] /= totalMass;
    linearMomentum[2] /= totalMass;

    if (print) {
      String sb = format(
          "  Center of Mass   (%12.3f,%12.3f,%12.3f)\n  Linear Momentum  (%12.3f,%12.3f,%12.3f)\n  Angular Momentum (%12.3f,%12.3f,%12.3f)",
          centerOfMass[0], centerOfMass[1], centerOfMass[2], linearMomentum[0], linearMomentum[1],
          linearMomentum[2], angularMomentum[0], angularMomentum[1], angularMomentum[2]);
      logger.info(sb);
    }

    if (remove) {
      removeCenterOfMassMotion(print);
      centerOfMassMotion(false, print);
    }
  }

  /**
   * Compute the current temperature and kinetic energy of the system.
   */
  public final void computeKineticEnergy() {
    double e = 0.0;
    double[] v = state.v();
    double[] mass = state.getMass();
    for (int i = 0; i < state.getNumberOfVariables(); i++) {
      double m = mass[i];
      if (m > 0.0) {
        double velocity = v[i];
        double v2 = velocity * velocity;
        e += m * v2;
      }
    }
    state.setTemperature(e / (kB * degreesOfFreedom));
    e *= 0.5 / KCAL_TO_GRAM_ANG2_PER_PS2;
    state.setKineticEnergy(e);
  }

  /**
   * The half-step temperature correction.
   *
   * @param dt a double.
   */
  public abstract void halfStep(double dt);

  /**
   * The full-step temperature correction.
   *
   * @param dt a double.
   */
  public abstract void fullStep(double dt);

  /**
   * Get the current temperature.
   *
   * <p>This depends on a previous call to the computeKineticEnergy.
   *
   * @return Current temperature.
   */
  public double getCurrentTemperature() {
    return state.getTemperature();
  }

  /**
   * Return the number of degrees of freedom.
   *
   * @return Degrees of freedom.
   */
  public int getDegreesOfFreedom() {
    return degreesOfFreedom;
  }

  /**
   * Get the current kinetic energy.
   *
   * <p>This depends on a previous call to the computeKineticEnergy.
   *
   * @return Kinetic energy.
   */
  public double getKineticEnergy() {
    return state.getKineticEnergy();
  }

  /**
   * Getter for the field <code>removeCenterOfMassMotion</code>.
   *
   * @return a boolean.
   */
  public boolean getRemoveCenterOfMassMotion() {
    return removeCenterOfMassMotion;
  }

  /**
   * If center of mass motion is being removed, then the mean kinetic energy of the system will be 3
   * kT/2 less than if center of mass motion is allowed.
   *
   * @param remove <code>true</code> if center of mass motion is being removed.
   */
  public void setRemoveCenterOfMassMotion(boolean remove) {
    removeCenterOfMassMotion = remove;
    int nVariables = state.getNumberOfVariables();
    if (removeCenterOfMassMotion) {
      degreesOfFreedom = nVariables - 3 - constrainedDoF;
    } else {
      degreesOfFreedom = nVariables - constrainedDoF;
    }
  }

  /**
   * Get the target temperature.
   *
   * @return Target temperature.
   */
  public double getTargetTemperature() {
    return targetTemperature;
  }

  /**
   * Set the target temperature.
   *
   * @param t Target temperature must be greater than absolute zero.
   * @since 1.0
   */
  public void setTargetTemperature(double t) {
    // Obey the Third Law of Thermodynamics.
    assert (t > 0.0);
    targetTemperature = t;
    kT = t * kB;
  }

  /**
   * Reset velocities from a Maxwell-Boltzmann distribution of momenta based on the supplied target
   * temperature. The variance of each independent momentum component is kT * mass.
   *
   * @param targetTemperature the target Temperature for the Maxwell distribution.
   */
  public void maxwell(double targetTemperature) {
    if (logger.isLoggable(Level.FINE)) {
      logger.fine("\n Initializing velocities to target temperature");
    }

    setTargetTemperature(targetTemperature);

    double[] v = state.v();
    double[] mass = state.getMass();
    for (int i = 0; i < state.getNumberOfVariables(); i++) {
      double m = mass[i];
      if (m > 0.0) {
        v[i] = random.nextGaussian() * sqrt(kB * targetTemperature / m);
      }
    }

    // Remove the center of mass motion.
    if (removeCenterOfMassMotion) {
      centerOfMassMotion(true, !quiet);
    }

    // Find the current kinetic energy and temperature.
    computeKineticEnergy();

    /*
     The current temperature will deviate slightly from the target
     temperature if the center of mass motion was removed and/or due to
     finite system size.

     Scale the velocities to enforce the target temperature.
    */
    double scale = sqrt(targetTemperature / state.getTemperature());
    for (int i = 0; i < state.getNumberOfVariables(); i++) {
      double m = mass[i];
      if (m > 0.0) {
        v[i] *= scale;
      }
    }

    // Update the kinetic energy and current temperature.
    computeKineticEnergy();

    log(Level.INFO);
  }

  /**
   * Reset velocities from a Maxwell-Boltzmann distribution based on the current target temperature
   * of thermostat.
   */
  public void maxwell() {
    maxwell(targetTemperature);
  }

  /**
   * Return 3 velocities from a Maxwell-Boltzmann distribution of momenta. The variance of each
   * independent momentum component is kT * mass.
   *
   * @param mass The mass for the degrees of freedom.
   * @return three velocity components.
   */
  public double[] maxwellIndividual(double mass) {
    double[] vv = new double[3];
    if (mass > 0.0) {
      for (int i = 0; i < 3; i++) {
        vv[i] = random.nextGaussian() * sqrt(kB * targetTemperature / mass);
      }
    }
    return vv;
  }

  /**
   * Setter for the field <code>quiet</code>.
   *
   * @param quiet a boolean.
   */
  public void setQuiet(boolean quiet) {
    this.quiet = quiet;
  }

  /**
   * The setRandomSeed method is used to initialize the Random number generator to the same starting
   * state, such that separate runs produce the same Maxwell-Boltzmann initial velocities. same
   *
   * @param seed The seed.
   */
  public void setRandomSeed(long seed) {
    random.setSeed(seed);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    return logTemp();
  }

  public String logTemp() {
    StringBuilder sb = new StringBuilder(
        format("  Target temperature:           %7.2f Kelvin\n", targetTemperature));
    sb.append(format("  Current temperature:          %7.2f Kelvin\n", state.getTemperature()));
    sb.append(format("  Number of variables:          %7d\n", state.getNumberOfVariables()));
    sb.append(format("  Number of degrees of freedom: %7d\n", degreesOfFreedom));
    sb.append(format("  Kinetic Energy:               %7.2f\n", state.getKineticEnergy()));
    sb.append(format("  kT per degree of freedom:     %7.2f",
        KCAL_TO_GRAM_ANG2_PER_PS2 * state.getKineticEnergy() / (degreesOfFreedom * kT)));
    return sb.toString();
  }

  /**
   * Log the target temperature and current number of kT per degree of freedom (should be 0.5 kT at
   * equilibrium).
   *
   * @param level a {@link java.util.logging.Level} object.
   */
  protected void log(Level level) {
    if (logger.isLoggable(level) && !quiet) {
      logger.log(level, logTemp());
    }
  }

  /**
   * Remove center of mass translational and rotational velocity by inverting the moment of inertia
   * tensor.
   *
   * @param print If true, log removal of center of mass motion.
   */
  private void removeCenterOfMassMotion(boolean print) {
    double xx = 0.0;
    double yy = 0.0;
    double zz = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yz = 0.0;
    int index = 0;
    double[] x = state.x();
    double[] mass = state.getMass();
    while (index < state.getNumberOfVariables()) {
      if (type[index] == VARIABLE_TYPE.OTHER) {
        index++;
        continue;
      }
      double m = mass[index];
      assert (type[index] == VARIABLE_TYPE.X);
      double xi = x[index++] - centerOfMass[0];
      assert (type[index] == VARIABLE_TYPE.Y);
      double yi = x[index++] - centerOfMass[1];
      assert (type[index] == VARIABLE_TYPE.Z);
      double zi = x[index++] - centerOfMass[2];
      xx += xi * xi * m;
      yy += yi * yi * m;
      zz += zi * zi * m;
      xy += xi * yi * m;
      xz += xi * zi * m;
      yz += yi * zi * m;
    }

    //        RealMatrix inertia = new Array2DRowRealMatrix(3, 3);
    //        inertia.setEntry(0, 0, yy + zz);
    //        inertia.setEntry(1, 0, -xy);
    //        inertia.setEntry(2, 0, -xz);
    //        inertia.setEntry(0, 1, -xy);
    //        inertia.setEntry(1, 1, xx + zz);
    //        inertia.setEntry(2, 1, -yz);
    //        inertia.setEntry(0, 2, -xz);
    //        inertia.setEntry(1, 2, -yz);
    //        inertia.setEntry(2, 2, xx + yy);
    //        inertia = new LUDecomposition(inertia).getSolver().getInverse();
    //        xx = inertia.getEntry(0, 0);
    //        yy = inertia.getEntry(1, 1);
    //        zz = inertia.getEntry(2, 2);
    //        xy = inertia.getEntry(0, 1);
    //        xz = inertia.getEntry(0, 2);
    //        yz = inertia.getEntry(1, 2);
    //        double ox = angularMomentum[0] * xx + angularMomentum[1] * xy + angularMomentum[2] *
    // xz;
    //        double oy = angularMomentum[0] * xy + angularMomentum[1] * yy + angularMomentum[2] *
    // yz;
    //        double oz = angularMomentum[0] * xz + angularMomentum[1] * yz + angularMomentum[2] *
    // zz;

    // Remove center of mass translational momentum.
    index = 0;
    double[] v = state.v();
    while (index < state.getNumberOfVariables()) {
      if (type[index] == VARIABLE_TYPE.OTHER) {
        index++;
        continue;
      }
      double m = mass[index];
      if (m > 0.0) {
        v[index++] -= linearMomentum[0] / m;
        v[index++] -= linearMomentum[1] / m;
        v[index++] -= linearMomentum[2] / m;
      } else {
        index += 3;
      }
    }

    // Only remove center of mass rotational momentum for non-periodic systems.
    //        if (false) {
    //            index = 0;
    //            while (index < nVariables) {
    //                if (type[index] == VARIABLE_TYPE.OTHER) {
    //                    index++;
    //                    continue;
    //                }
    //                double xi = x[index++] - centerOfMass[0];
    //                double yi = x[index++] - centerOfMass[1];
    //                double zi = x[index] - centerOfMass[2];
    //                index -= 2;
    //                v[index++] += (-oy * zi + oz * yi);
    //                v[index++] += (-oz * xi + ox * zi);
    //                v[index++] += (-ox * yi + oy * xi);
    //            }
    //        }

    if (print) {
      logger.info("  Center of mass motion removed.");
    }
  }
}
