//******************************************************************************
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
//******************************************************************************
package ffx.potential.groovy

import edu.rit.pj.ParallelTeam
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import ffx.potential.nonbonded.implicit.ConnollyRegion
import ffx.potential.nonbonded.implicit.GaussVol
import ffx.potential.parameters.ForceField
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static ffx.utilities.Constants.NS2SEC
import static java.lang.String.format
import static org.apache.commons.math3.util.FastMath.pow

/**
 * Calculate the surface area and volume using the GaussVol (default) or Connolly algorithm.
 * <br>
 * Usage:
 * <br>
 * ffxc Volume [options] &lt;filename&gt;
 */
@Command(description = " Calculate the surface area and volume using the GaussVol (default) or Connolly algorithm.",
    name = "Volume")
class Volume extends PotentialScript {

  private static final double rminToSigma = 1.0 / pow((double) 2.0, (double) 1.0 / 6.0)

  /**
   * -c or --connolly Use the Connolly algorithm to compute volume and surface area (instead of GaussVol).
   */
  @Option(names = ['-c', '--connolly'], paramLabel = "false",
      description = "Use the Connolly algorithm to compute solvent excluded volume and solvent accessible surface area.")
  private boolean connolly = false

  /**
   * -m or --molecular For Connolly, compute molecular volume and surface area (instead of SEV/SASA).
   */
  @Option(names = ['-m', '--molecular'], paramLabel = "false",
      description = "For Connolly, compute molecular volume and surface area (instead of SEV/SASA).")
  private boolean molecular = false

  /**
   * --vdW or --vanDerWaals For Connolly, compute van der Waals volume and surface area (instead of SEV/SASA).
   */
  @Option(names = ['--vdW', '--vanDerWaals'], paramLabel = "false",
      description = "For Connolly, compute van der Waals volume and surface area (instead of SEV/SASA)")
  private boolean vdW = false

  /**
   * -p or --probe For Connolly, set the exclude radius (SASA) or probe (molecular surface). Ignored for vdW.
   */
  @Option(names = ['-p', '--probe'], paramLabel = "1.4",
      description = "For Connolly, set the exclude radius (SASA) or probe radius (molecular surface). Ignored for vdW.")
  private double probe = 1.4

  /**
   * -y or --includeHydrogen Include Hydrogen in calculation volume and surface area.
   */
  @Option(names = ['-y', '--includeHydrogen'], paramLabel = "false",
      description = "Include Hydrogen in calculation volume and surface area.")
  private boolean includeHydrogen = false

  /**
   * -s or --sigma Use sigma radii instead of Rmin.
   */
  @Option(names = ['-s', '--sigma'], paramLabel = "false",
      description = "Use sigma radii instead of Rmin.")
  private boolean sigma = false

  /**
   * -o or --offset For GaussVol, add an offset to all atomic radii.
   */
  @Option(names = ['-o', '--offset'], paramLabel = "0.0",
      description = "Add an offset to all atomic radii for GaussVol volume and surface area.")
  private double offset = 0.0

  /**
   * -v or --verbose enables printing out all energy components for multi-snapshot files (
   * the first snapshot is always printed verbosely).
   */
  @Option(names = ['-v', '--verbose'], paramLabel = "false",
      description = "Print out all components of volume of molecule and offset.")
  private boolean verbose = false


  /**
   * The final argument is an atomic coordinate file in PDB or XYZ format.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  String filename = null

  /**
   * JUnit Testing Variables
   */
  public double totalVolume = 0.0
  public double totalSurfaceArea = 0.0

  /**
   * Volume Constructor.
   */
  Volume() {
    this(new Binding())
  }

  /**
   * Volume Constructor.
   * @param binding Groovy Binding to use.
   */
  Volume(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  Volume run() {

    // Init the context and bind variables.
    if (!init()) {
      return null
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return null
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    logger.info("\n Calculating volume and surface area for " + filename)

    Atom[] atoms = activeAssembly.getAtomArray()
    int nAtoms = atoms.length

    if (!connolly) {
      // Input
      double[][] positions = new double[nAtoms][3]
      int index = 0
      for (Atom atom : atoms) {
        positions[index][0] = atom.getX()
        positions[index][1] = atom.getY()
        positions[index][2] = atom.getZ()
        index++
      }

      ForceField forceField = activeAssembly.getForceField()
      if (includeHydrogen) {
        forceField.addProperty("GAUSSVOL_HYDROGEN", "true")
      }
      if (sigma) {
        forceField.addProperty("GAUSSVOL_USE_SIGMA", "true")
      }
      if (offset > 0.0) {
        forceField.addProperty("GAUSSVOL_RADII_OFFSET", Double.toString(offset))
      }

      if (!forceField.hasProperty("GAUSSVOL_RADII_SCALE")) {
        forceField.addProperty("GAUSSVOL_RADII_SCALE", Double.toString(1.0))
      }


      // Run Volume calculation to get vdw volume of molecule
      ParallelTeam parallelTeam = new ParallelTeam()
      GaussVol gaussVol = new GaussVol(atoms, activeAssembly.getForceField(), parallelTeam)
      long time = -System.nanoTime()
      gaussVol.computeVolumeAndSA(positions)
      time += System.nanoTime()

      if (verbose) {
        logger.info(format("\n Maximum depth of overlaps in tree: %d", gaussVol.getMaximumDepth()))
        logger.info(format(" Total number of overlaps in tree: %d", gaussVol.getTotalNumberOfOverlaps()))

        double[] radii = gaussVol.getRadii()
        index = 0
        for (Atom atom : atoms) {
          logger.info(format(" %5d %4s: %6.4f", index, atom.getName(), radii[index]))
          index++
        }
      }

      logger.info("\n GaussVol Surface Area and Volume\n")
      if (sigma) {
        logger.info(format("  Radii:                  Sigma"))
      } else {
        logger.info(format("  Radii:                   Rmin"))
      }
      logger.info(format("  Radii offset:        %8.4f (Ang)", offset))
      logger.info(format("  Include hydrogen:    %8b", includeHydrogen))
      logger.info(format("  Volume:            %10.4f (Ang^3)", gaussVol.getVolume()))
      logger.info(format("  Surface Area:      %10.4f (Ang^2)", gaussVol.getSurfaceArea()))
      logger.info(format("  Time:              %10.4f (sec)", time * NS2SEC))

      // Set JUnit testing variables based on output volume and surface area
      totalVolume = gaussVol.getVolume()
      totalSurfaceArea = gaussVol.getSurfaceArea()
    } else {
      // For Connolly molecular volume & surface area, use the chosen probe and set exclude to 0.0.
      double exclude = 0.0

      if (vdW) {
        // For Connolly vdW, both exclude & probe are zero.
        exclude = 0.0
        probe = 0.0
      } else if (!molecular) {
        // For Connolly SEV/SASA, set exclude to the chosen probe, and zero the probe.
        exclude = probe
        probe = 0.0
      }

      double[] radii = new double[nAtoms]
      int index = 0
      for (Atom atom : atoms) {
        radii[index] = atom.getVDWType().radius / 2.0
        if (sigma) {
          radii[index] *= rminToSigma
        }
        boolean hydrogen = atom.isHydrogen()
        if (!includeHydrogen && hydrogen) {
          radii[index] = 0.0
        }
        index++
      }

      // Note that the VolumeRegion code is currently limited to a single thread.
      int nThreads = 1
      ConnollyRegion connollyRegion = new ConnollyRegion(atoms, radii, nThreads)
      // For solvent excluded volume.
      connollyRegion.setExclude(exclude)
      // For molecular surface.
      connollyRegion.setProbe(probe)
      long time = -System.nanoTime()
      connollyRegion.runVolume()
      time += System.nanoTime()

      double zero = 0.0
      if (vdW || (probe == zero && exclude == zero)) {
        logger.info("\n Connolly van der Waals Surface Area and Volume\n")
      } else if (!molecular) {
        logger.info("\n Connolly Solvent Accessible Surface Area and Solvent Excluded Volume\n")
        logger.info(format("  Exclude radius:      %8.4f (Ang)", exclude))
      } else {
        logger.info("\n Connolly Molecular Surface Area and Volume")
        logger.info(format("  Probe radius:       %8.4f (Ang)", probe))
      }
      if (sigma) {
        logger.info(format("  Radii:                  Sigma"))
      } else {
        logger.info(format("  Radii:                   Rmin"))
      }
      logger.info(format("  Include hydrogen:    %8b", includeHydrogen))
      logger.info(format("  Volume:            %10.4f (Ang^3)", connollyRegion.getVolume()))
      logger.info(format("  Surface Area:      %10.4f (Ang^2)", connollyRegion.getSurfaceArea()))
      logger.info(format("  Time:              %10.4f (sec)", time * NS2SEC))

      // Set JUnit testing variables based on output volume and surface area
      totalVolume = connollyRegion.getVolume()
      totalSurfaceArea = connollyRegion.getSurfaceArea()
    }

    return this
  }
}
