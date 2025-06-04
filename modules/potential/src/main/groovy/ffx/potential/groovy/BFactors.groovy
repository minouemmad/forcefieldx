
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

import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.potential.cli.PotentialScript
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters
import picocli.CommandLine.Option
import org.apache.commons.io.FilenameUtils

@Command(description = " Extract B-factors/confidence scores for ALL atoms from a PDB file", name = "BFactors")
class BFactors extends PotentialScript {

    @Option(names = ["-c", "--csv"], paramLabel = "output.csv",
            description = "Optional CSV file to append B-factors to. If not provided, a new CSV will be created.")
    private String csvFilename = null

    @Option(names = ["-d", "--delimiter"], paramLabel = ",",
            description = "Delimiter for CSV file (default is comma)")
    private String delimiter = ","

    @Parameters(arity = "1", paramLabel = "file",
            description = 'The atomic coordinate file in PDB format')
    private String filename = null

    BFactors() { this(new Binding()) }
    BFactors(Binding binding) { super(binding) }

    @Override
    BFactors run() {
        if (!init()) return null
        activeAssembly = getActiveAssembly(filename)
        if (activeAssembly == null) {
            logger.info(helpString())
            return null
        }

        List<Residue> residues = activeAssembly.getResidueList()
        List<Map<String, Object>> bFactorData = collectBFactors(residues)
        writeBFactorsToCSV(bFactorData)

        // Calculate and print median B-factor
        printMedianBFactor(bFactorData)

        return this
    }

    private List<Map<String, Object>> collectBFactors(List<Residue> residues) {
        List<Map<String, Object>> bFactorData = []

        residues.each { residue ->
            residue.getAtomList().each { atom ->
                bFactorData << [
                        resNum: residue.getResidueNumber(),
                        resName: residue.getName(),
                        atomName: atom.getName(),
                        chain: residue.getChainID(),
                        bFactor: atom.getTempFactor()
                ]
            }
        }
        return bFactorData
    }

    private void writeBFactorsToCSV(List<Map<String, Object>> bFactorData) {
        String baseName = FilenameUtils.getBaseName(filename)
        String outputFilename = csvFilename ?: "b_factors_${baseName}.csv"

        try {
            File outputFile = new File(outputFilename)
            boolean fileExists = outputFile.exists()

            outputFile.withWriter { writer ->
                if (!fileExists || !csvFilename) {
                    writer.write("Residue Number${delimiter}Residue Name${delimiter}Atom Name${delimiter}Chain${delimiter}B-factor\n")
                }

                bFactorData.each { data ->
                    writer.write("${data.resNum}${delimiter}${data.resName}${delimiter}${data.atomName}${delimiter}${data.chain}${delimiter}${data.bFactor}\n")
                }
            }

            logger.info("\n Successfully wrote B-factors to ${outputFilename}")
        } catch (Exception e) {
            logger.severe("Error writing B-factors to CSV: ${e.message}")
        }
    }

    private void printMedianBFactor(List<Map<String, Object>> bFactorData) {
        // Extract all B-factors
        List<Double> bFactors = bFactorData.collect { it.bFactor }

        // Sort the list
        bFactors.sort()

        // Calculate median
        double median
        int size = bFactors.size()
        if (size == 0) {
            logger.info("\n No B-factors found to calculate median")
            return
        }

        if (size % 2 == 0) {
            median = (bFactors[size/2 - 1] + bFactors[size/2]) / 2.0
        } else {
            median = bFactors[size/2]
        }

        // Print results
        logger.info("\n B-factor statistics:")
        logger.info(String.format("  Total atoms: %d", size))
        logger.info(String.format("  Minimum B-factor: %.2f", bFactors.min()))
        logger.info(String.format("  Maximum B-factor: %.2f", bFactors.max()))
        logger.info(String.format("  Median B-factor: %.2f", median))
    }
}