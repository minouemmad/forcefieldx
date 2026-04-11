// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.potential.commands;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.cli.PotentialCommand;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.utilities.FFXBinding;
import picocli.CommandLine.Command;
import picocli.CommandLine.Parameters;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import static java.lang.String.format;

/**
 * Evaluate the AMOEBA energy of a system and write the permanent multipoles
 * (charge, dipole x/y/z, quadrupole xx/xy/yy/xz/yz/zz) for every atom to a
 * {@code .uperm} file in the same directory as the input structure.
 *
 * <p>The 10 components of {@code globalMultipole[0][atom]} follow the FFX/Tinker convention:
 * <pre>
 *   [0]  charge (monopole)
 *   [1]  dipole x
 *   [2]  dipole y
 *   [3]  dipole z
 *   [4]  quadrupole xx
 *   [5]  quadrupole xy
 *   [6]  quadrupole yy
 *   [7]  quadrupole xz
 *   [8]  quadrupole yz
 *   [9]  quadrupole zz
 * </pre>
 *
 * <p>Usage:
 * <pre>
 *   ffxc SavePermanentMoments &lt;filename.pdb&gt;
 * </pre>
 */
@Command(name = "SavePermanentMoments",
    description = " Evaluate AMOEBA energy and save permanent multipoles to a .uperm file.")
public class SavePermanentMoments extends PotentialCommand {

  @Parameters(arity = "1", paramLabel = "file",
      description = "PDB or XYZ input file.")
  private String filename = null;

  public SavePermanentMoments() {
    super();
  }

  public SavePermanentMoments(FFXBinding binding) {
    super(binding);
  }

  public SavePermanentMoments(String[] args) {
    super(args);
  }

  @Override
  public SavePermanentMoments run() {
    if (!init()) {
      return this;
    }

    activeAssembly = getActiveAssembly(filename);
    if (activeAssembly == null) {
      logger.info(helpString());
      return this;
    }

    filename = activeAssembly.getFile().getAbsolutePath();
    logger.info(format("\n Running SavePermanentMoments on %s", filename));

    // Evaluate the energy to populate globalMultipole.
    ForceFieldEnergy energy = activeAssembly.getPotentialEnergy();
    double[] x = new double[energy.getNumberOfVariables()];
    energy.getCoordinates(x);
    energy.energy(x, true);

    // Retrieve the PME node.
    ParticleMeshEwald pme = energy.getPmeNode();
    if (pme == null || pme.globalMultipole == null || pme.globalMultipole.length == 0) {
      logger.warning(" No PME permanent multipoles available — is an AMOEBA force field active?");
      return this;
    }

    double[][] multipoles = pme.globalMultipole[0];
    Atom[] atoms = activeAssembly.getAtomArray();

    if (multipoles == null || multipoles.length < atoms.length) {
      logger.warning(" Incomplete globalMultipole array; skipping output.");
      return this;
    }

    // Build output path: same directory as input, .uperm extension.
    File inputFile = new File(filename);
    String baseName = inputFile.getName();
    int dot = baseName.lastIndexOf('.');
    if (dot > 0) {
      baseName = baseName.substring(0, dot);
    }
    File outFile = new File(inputFile.getParentFile(), baseName + ".uperm");

    try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFile))) {
      writer.write(format("%7d  %s%n", atoms.length, activeAssembly.getName()));
      for (int i = 0; i < atoms.length; i++) {
        Atom a = atoms[i];
        double[] m = multipoles[i];
        // charge  dipX dipY dipZ  quadXX quadXY quadYY quadXZ quadYZ quadZZ
        writer.write(format(
            "%7d %-4s %14.8f  %14.8f %14.8f %14.8f  %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f%n",
            a.getIndex(), a.getName(),
            m[0],
            m[1], m[2], m[3],
            m[4], m[5], m[6], m[7], m[8], m[9]));
      }
    } catch (IOException ex) {
      logger.warning(format(" Exception writing permanent multipoles to %s: %s", outFile, ex));
      return this;
    }

    logger.info(format(" Wrote permanent multipoles to: %s", outFile.getAbsolutePath()));
    return this;
  }

  @Override
  public List<Potential> getPotentials() {
    if (activeAssembly == null) {
      return Collections.emptyList();
    }
    ForceFieldEnergy ffe = activeAssembly.getPotentialEnergy();
    return ffe != null ? Collections.singletonList(ffe) : Collections.emptyList();
  }
}
