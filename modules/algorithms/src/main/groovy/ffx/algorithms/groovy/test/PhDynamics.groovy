//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.algorithms.groovy.test

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.dynamics.MolecularDynamics
import ffx.algorithms.ph.PhMD
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Residue
import ffx.potential.cli.WriteoutOptions
import ffx.potential.extended.TitrationUtils
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The Dynamics script implements molecular and stochastic dynamics algorithms.
 * <br>
 * Usage:
 * <br>
 * ffxc CPHMDynamics [options] &lt;filename&gt; [file2...]
 */
@Command(description = " Run constant pH dynamics on a system.", name = "ffxc PHDynamics")
class PhDynamics extends AlgorithmsScript {

  @Mixin
  DynamicsOptions dynamics

  @Mixin
  WriteoutOptions writeOut

  /**
   * --pH or --constantPH Constant pH value for molecular dynamics.
   */
  @Option(names = ['--pH', '--constantPH'], paramLabel = '7.4',
      description = 'Constant pH value for molecular dynamics')
  double pH = 7.4

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = "XYZ or PDB input files.")
  private List<String> filenames

  /**
   * Creation of a public field to try and make the JUnit test work, original code does not declare this as a public field.
   * Originally it is declared in the run method
   */
  public Potential potential = null
  public MolecularDynamics molDyn = null

  MolecularDynamics getMolecularDynamics() {
    return molDyn
  }

  Potential getPotentialObject() {
    return potential
  }

  /**
   * Dynamics Constructor.
   */
  PhDynamics() {
    this(new Binding())
  }

  /**
   * Dynamics Constructor.
   * @param binding The Groovy Binding to use.
   */
  PhDynamics(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  PhDynamics run() {

    if (!init()) {
      return this
    }

    dynamics.init()

    TitrationUtils.initDiscountPreloadProperties()

    String modelFilename
    if (filenames != null && filenames.size() > 0) {
      MolecularAssembly molecularAssembly =
          TitrationUtils.openFullyProtonated(new File(filenames.get(0)))
      MolecularAssembly[] assemblies = [molecularAssembly]
      activeAssembly = assemblies[0]
      modelFilename = filenames.get(0)
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      modelFilename = activeAssembly.getFile().getAbsolutePath()
    }

    // Select all possible titrating residues.
    List<Residue> titrating = TitrationUtils.chooseTitratables(activeAssembly)

    potential = activeAssembly.getPotentialEnergy()
    double[] x = new double[potential.getNumberOfVariables()]
    potential.getCoordinates(x)
    potential.energy(x, true)

    logger.info("\n Running molecular dynamics on " + modelFilename)

    molDyn = dynamics.getDynamics(writeOut, potential, activeAssembly, algorithmListener)

    PhMD.Distribution distribution = PhMD.Distribution.CONTINUOUS
    // No MC steps yet.
    int mcStepFrequency = 10
    PhMD phmd = new PhMD(distribution, activeAssembly, molDyn, titrating, pH, mcStepFrequency)

    // Restart File
    File dyn = new File(FilenameUtils.removeExtension(modelFilename) + ".dyn")
    if (!dyn.exists()) {
      dyn = null
    }

    molDyn.dynamic(dynamics.steps, dynamics.dt,
        dynamics.report, dynamics.write, dynamics.temperature, true, dyn)

    return this
  }

  /**
   * {@inheritDoc}
   */
  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (potential == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(potential)
    }
    return potentials
  }

}