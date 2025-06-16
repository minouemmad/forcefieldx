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
import ffx.potential.utils.GetProteinFeatures
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters
import picocli.CommandLine.Option
import org.apache.commons.io.FilenameUtils

@Command(description = " Extract B-factors/confidence scores and structural features from PDB files", name = "BFactors")
class BFactors extends PotentialScript {

    @Option(names = ["-c", "--csv"], paramLabel = "output.csv",
            description = "Optional CSV file to append B-factors to. If not provided, a new CSV will be created.")
    private String csvFilename = null

    @Option(names = ["-d", "--delimiter"], paramLabel = ",",
            description = "Delimiter for CSV file (default is comma)")
    private String delimiter = ","

    @Option(names = ["-b", "--batch"], paramLabel = "false",
            description = "Process a directory of PDB files and output features to PDB_Features.csv")
    private boolean batchMode = false

    @Option(names = ["-m", "--maestro"], paramLabel = "false",
            description = "Calculate Instability Heat Score (IHS) using MAESTRO predictions")
    private boolean calculateIHS = false

    @Option(names = ["-s", "--stabilityDir"], paramLabel = "stabilityPredictions",
            description = "Directory containing MAESTRO prediction files (.maestro.txt)")
    private String stabilityDir = "stabilityPredictions"

    @Parameters(arity = "1", paramLabel = "file",
            description = 'The atomic coordinate file in PDB format or directory containing PDB files')
    private String filename = null

    BFactors() { this(new Binding()) }
    BFactors(Binding binding) { super(binding) }

    @Override
    BFactors run() {
        if (!init()) return null

        if (batchMode) {
            processDirectory()
        } else {
            processSingleFile()
        }

        return this
    }

    private void processSingleFile() {
        activeAssembly = getActiveAssembly(filename)
        if (activeAssembly == null) {
            logger.info(helpString())
            return
        }

        List<Residue> residues = activeAssembly.getResidueList()
        List<Map<String, Object>> bFactorData = collectBFactors(residues)
        writeBFactorsToCSV(bFactorData)

        // Calculate and print median B-factor
        printMedianBFactor(bFactorData)

        // Add IHS calculation for single file
        if (calculateIHS) {
            String pdbId = FilenameUtils.getBaseName(filename)
            double ihs = calculateIHS(residues, pdbId)
            logger.info("\nFinal IHS for ${pdbId}: ${String.format("%.2f", ihs)}\n")

            // You could also write this to a CSV if needed
            Map<String, Object> result = [
                    pdbId: pdbId,
                    medianBfactor: String.format("%.2f", calculateMedian(bFactorData)),
                    totalSASA: "0.0", // Would need to calculate SASA as done in processDirectory()
                    numResidues: residues.size(),
                    IHS: String.format("%.2f", ihs)
            ]
            writeExtendedResultsToCSV([result])
        }
    }

    private void processDirectory() {
        File directory = new File(filename)
        if (!directory.isDirectory()) {
            logger.severe("The provided path is not a directory: ${filename}")
            return
        }

        List<Map<String, Object>> fileResults = []
        List<Map<String, Object>> geneHeatScores = []
        File[] pdbFiles = directory.listFiles({ File f -> f.name.toLowerCase().endsWith('.pdb') } as FileFilter)

        if (pdbFiles == null || pdbFiles.length == 0) {
            logger.severe("No PDB files found in directory: ${filename}")
            return
        }

        // Create PDB to gene mapping
        Map<String, String> pdbToGene = createPdbToGeneMap()

        logger.info("\nProcessing ${pdbFiles.length} PDB files in directory: ${filename}")
        int processedCount = 0
        int errorCount = 0

        // Check stability directory exists if IHS calculation is requested
        File stabilityDirectory = null
        if (calculateIHS) {
            stabilityDirectory = new File(stabilityDir)
            if (!stabilityDirectory.exists() || !stabilityDirectory.isDirectory()) {
                logger.severe("Stability predictions directory not found: ${stabilityDir}")
                return
            }
        }

        pdbFiles.each { pdbFile ->
            try {
                logger.info("Processing file: ${pdbFile.name}")

                // Create a new binding with required properties
                Binding fileBinding = new Binding()
                fileBinding.setVariable("args", [pdbFile.absolutePath] as String[])
                fileBinding.setVariable("csvFilename", csvFilename)
                fileBinding.setVariable("delimiter", delimiter)
                fileBinding.setVariable("batchMode", batchMode)
                fileBinding.setVariable("calculateIHS", calculateIHS)
                fileBinding.setVariable("stabilityDir", stabilityDir)

                // Create new script instance with initialized binding
                BFactors fileProcessor = new BFactors(fileBinding)
                fileProcessor.run()

                // Check if assembly was created successfully
                if (fileProcessor.activeAssembly == null) {
                    logger.warning("Could not create assembly for file: ${pdbFile.name}")
                    errorCount++
                    return // continues to next file
                }

                // Get structural features
                List<Residue> residues = fileProcessor.activeAssembly.residueList
                List<Map<String, Object>> bFactorData = fileProcessor.collectBFactors(residues)

                if (bFactorData.isEmpty()) {
                    logger.warning("No B-factors collected from file: ${pdbFile.name}")
                    errorCount++
                    return
                }

                // Initialize GetProteinFeatures and calculate SASA
                GetProteinFeatures proteinFeatures = new GetProteinFeatures()
                residues.each { residue ->
                    double residueSA = GetProteinFeatures.standardSurfaceArea.getOrDefault(
                            residue.getAminoAcid3(), 0.0)
                    proteinFeatures.saveFeatures(residue, residueSA, false, false)
                }
                double totalSASA = proteinFeatures.getTotalSurfaceArea()

                // Calculate median B-factor
                double median = calculateMedian(bFactorData)
                String pdbId = FilenameUtils.getBaseName(pdbFile.name)

                // Calculate IHS if requested
                double ihs = 0.0
                if (calculateIHS) {
                    ihs = calculateIHS(residues, pdbId)
                    logger.info("\nFinal IHS for ${pdbId}: ${String.format("%.2f", ihs)}\n")
                }

                // Collect all results
                Map<String, Object> result = [
                        pdbId: pdbId,
                        medianBfactor: String.format("%.2f", median),
                        totalSASA: String.format("%.2f", totalSASA),
                        numResidues: residues.size()
                ]

                if (calculateIHS) {
                    result["IHS"] = String.format("%.2f", ihs)
                }

                // After calculating IHS, add to geneHeatScores if available
                if (calculateIHS) {
                    String cleanPdbId = FilenameUtils.getBaseName(pdbFile.name).toUpperCase()
                    String geneName = pdbToGene.get(cleanPdbId, "Unknown")

                    geneHeatScores << [
                            pdbId: cleanPdbId,
                            geneName: geneName,
                            heatScore: String.format("%.2f", ihs)
                    ]
                }

                fileResults << result

                logger.info(String.format("  Processed %s - Median B-factor: %.2f, Total SASA: %.2f, Residues: %d",
                        pdbId, median, totalSASA, residues.size()))
                processedCount++

            } catch (Exception e) {
                logger.warning("Error processing file ${pdbFile.name}: ${e.message}")
                errorCount++
            }
        }

        writeExtendedResultsToCSV(fileResults)
        // Write gene heat scores to separate CSV if IHS was calculated
        if (calculateIHS && !geneHeatScores.isEmpty()) {
            writeGeneHeatScoresToCSV(geneHeatScores)
        }

        logger.info("\nProcessing complete:")
        logger.info(String.format("  Successfully processed: %d files", processedCount))
        logger.info(String.format("  Files with errors: %d", errorCount))
        logger.info("  Results written to PDB_Features.csv")
        if (calculateIHS) {
            logger.info("  Gene heat scores written to Gene_Heat_Scores.csv")
        }
    }

    private void writeGeneHeatScoresToCSV(List<Map<String, Object>> geneHeatScores) {
        String outputFilename = "Gene_Heat_Scores.csv"

        try {
            File outputFile = new File(outputFilename)
            outputFile.withWriter { writer ->
                // Write header
                writer.write("PDB_ID${delimiter}Gene_Name${delimiter}Heat_Score\n")

                // Write data sorted by heat score (descending)
                geneHeatScores.sort { -it.heatScore.toDouble() }.each { result ->
                    writer.write("${result.pdbId}${delimiter}${result.geneName}${delimiter}${result.heatScore}\n")
                }
            }
            logger.info("\nSuccessfully wrote gene heat scores to ${outputFilename}")
        } catch (Exception e) {
            logger.severe("Error writing gene heat scores to CSV: ${e.message}")
        }
    }
    private Map<String, String> createPdbToGeneMap() {
        return [
                "1WKU": "ACTN3",
                "2R0O": "ACTN4",
                "3MY0": "ACVRL1",
                "3PH9": "AGR3",
                "2BP1": "AKR7A2",
                "2X18": "AKT3",
                "1CPU": "AMY2A",
                "3L81": "AP4M1",
                "1E3G": "AR",
                "5NP0": "ATM",
                "2HLQ": "BMPR2",
                "1K2P": "BTK",
                "5C9J": "CD1C",
                "4AAA": "CDKL2",
                "1KWM": "CPB1",
                "2JDG": "CRYGB",
                "1YB5": "CRYZ",
                "1FH0": "CTSV",
                "2B5L": "CUL4B",
                "2I4I": "DDX3X",
                "3ZFM": "EPHB2",
                "5HYN": "EZH2",
                "3ZGQ": "FIT5",
                "2C6Y": "FOXK2",
                "2A2C": "GALK2",
                "1F5N": "GBP1",
                "1V4S": "GCK",
                "1A22": "GH1",
                "2QKH": "GIPR",
                "5ER7": "GJB2",
                "2H8R": "HNF1B",
                "1T0P": "ICAM3",
                "2M1W": "IL18RAP",
                "4UXT": "KIF5A",
                "3KXZ": "LCK",
                "1N7D": "LDLR",
                "2P4E": "LRP8",
                "5TJA": "MCOLN1",
                "3U84": "MEN1",
                "4XI6": "MIB1",
                "1JAP": "MMP8",
                "3THX": "MSH3",
                "4PA0": "MYH7",
                "1H4R": "NF2",
                "1IKN": "NFKBIA",
                "5UHV": "NRAS",
                "4FZV": "NSUN2",
                "3Q91": "NUDT14",
                "4NHX": "OGFOD1",
                "2DEW": "PADI4",
                "1DMW": "PAH",
                "1PDV": "PARK7",
                "3ELO": "PLA3G1B",
                "3IKM": "POLG",
                "4OR9": "PPP3CB",
                "5BZZ": "PTEN",
                "2SHP": "PTPN11",
                "4PJU": "RAD21",
                "1XAP": "RARB",
                "2F8X": "RBPJ",
                "4MTH": "REG3A",
                "6GXZ": "RPAP3",
                "4D9T": "RPS6KA3",
                "4UYB": "SEC14L3",
                "1ANT": "SERPINC1",
                "4MHX": "SGSH",
                "4PYP": "SLC2A1",
                "3VFD": "SPAST",
                "1YVL": "STAT1",
                "3BFX": "SULT1C2",
                "4LG9": "TBL1XR1",
                "3KAS": "TFRC",
                "2VFJ": "TNFAIP3",
                "4YFI": "TNNI3K",
                "1QBK": "TNPO1",
                "1A36": "TOP1",
                "1UOL": "TP53",
                "1FLK": "TRAF3",
                "4OLI": "TYK2",
                "1UOU": "TYMP",
                "1Y8Q": "UBA2",
                "3HU3": "VCP",
                "4NS5": "ZMYND11"
        ]
    }

    /**
     * Map PDB ID to UniProt ID using UniProt's REST API
     */
    private List<String> mapPdbToUniprot(String pdbId) {
        logger.info("Attempting to map PDB ${pdbId} to UniProt IDs using REST API...")
        List<String> uniprotIds = []

        try {
            // Step 1: Submit ID mapping job
            def submitUrl = "https://rest.uniprot.org/idmapping/run"
            def submitCmd = ["curl", "--request", "POST", submitUrl,
                             "--form", "ids=${pdbId}",
                             "--form", "from=PDB",
                             "--form", "to=UniProtKB"]

            logger.info("Executing: ${submitCmd.join(' ')}")
            def submitProcess = submitCmd.execute()
            submitProcess.waitFor()

            if (submitProcess.exitValue() != 0) {
                logger.warning("Submit command failed with exit code ${submitProcess.exitValue()}")
                logger.warning("Error output: ${submitProcess.err.text}")
                return uniprotIds
            }

            def submitResponse = submitProcess.text
            logger.info("Submit response: ${submitResponse}")

            // Parse job ID from response
            def jobId = (submitResponse =~ /"jobId":"([^"]+)"/)?.with { it[0][1] }
            if (!jobId) {
                logger.warning("Could not extract job ID from response")
                return uniprotIds
            }
            logger.info("Job ID: ${jobId}")

            // Step 2: Poll for job completion
            def statusUrl = "https://rest.uniprot.org/idmapping/status/${jobId}"
            int attempts = 0
            int maxAttempts = 10
            boolean jobComplete = false

            while (attempts < maxAttempts && !jobComplete) {
                sleep(3000) // Wait 3 seconds between checks
                attempts++

                def statusCmd = ["curl", "-s", statusUrl]
                def statusProcess = statusCmd.execute()
                statusProcess.waitFor()

                if (statusProcess.exitValue() != 0) {
                    logger.warning("Status check failed for attempt ${attempts}")
                    continue
                }

                def statusResponse = statusProcess.text
                if (statusResponse.contains('"jobStatus":"FINISHED"')) {
                    jobComplete = true
                }
            }

            if (!jobComplete) {
                logger.warning("Job did not complete within ${maxAttempts} attempts")
                return uniprotIds
            }

            // Step 3: Get results
            def resultsUrl = "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/${jobId}?format=tsv"
            def resultsCmd = ["curl", "-s", resultsUrl]
            logger.info("Executing: ${resultsCmd.join(' ')}")

            def resultsProcess = resultsCmd.execute()
            resultsProcess.waitFor()

            if (resultsProcess.exitValue() != 0) {
                logger.warning("Results command failed with exit code ${resultsProcess.exitValue()}")
                logger.warning("Error output: ${resultsProcess.err.text}")
                return uniprotIds
            }

            def results = resultsProcess.text
            logger.info("Results:\n${results}")

            // Parse TSV results - collect all UniProt IDs
            def lines = results.readLines()
            if (lines.size() > 1) {
                lines[1..-1].each { line ->
                    def parts = line.split("\t")
                    if (parts.size() >= 2) {
                        String uniprotId = parts[1].trim()
                        if (uniprotId && !uniprotIds.contains(uniprotId)) {
                            uniprotIds.add(uniprotId)
                        }
                    }
                }
            }

            if (uniprotIds) {
                logger.info("Successfully mapped to UniProt IDs: ${uniprotIds.join(', ')}")
            } else {
                logger.warning("No UniProt IDs found in results")
            }

        } catch (Exception e) {
            logger.warning("Error mapping PDB to UniProt: ${e.message}")
            e.printStackTrace()
        }

        return uniprotIds
    }

    /**
     * Find MAESTRO file using UniProt ID pattern with detailed logging
     */
    private File findMaestroFile(String pdbId, File maestroDir) {
        logger.info("\nSearching for MAESTRO file matching PDB: ${pdbId}")

        // First try direct PDB ID match
        logger.info("Attempting direct PDB ID match...")
        File directMatch = maestroDir.listFiles().find {
            it.name.startsWith(pdbId) && it.name.endsWith(".maestro.txt")
        }
        if (directMatch) {
            logger.info("Found direct match: ${directMatch.name}")
            return directMatch
        }

        // If no direct match, try UniProt mapping
        logger.info("No direct match found, attempting UniProt mapping...")
        List<String> uniprotIds = mapPdbToUniprot(pdbId)
        if (!uniprotIds) {
            logger.warning("Could not map PDB ID ${pdbId} to any UniProt IDs")
            return null
        }
        logger.info("Mapped ${pdbId} to UniProt IDs: ${uniprotIds.join(', ')}")

        // Try all possible patterns for each UniProt ID
        for (String uniprotId : uniprotIds) {
            logger.info("Searching for files matching UniProt ID: ${uniprotId}")

            // Define patterns for this specific uniprotId
            def patterns = [
                    ~/(?i)AF-${uniprotId}-F1.*\.maestro\.txt/,
                    ~/(?i)${uniprotId}\.maestro\.txt/,
                    ~/(?i)${uniprotId}-.*\.maestro\.txt/
            ]

            // Try all pattern variations
            for (pattern in patterns) {
                def potentialFiles = maestroDir.listFiles().findAll { file ->
                    file.name ==~ pattern
                }

                if (potentialFiles) {
                    logger.info("Found ${potentialFiles.size()} matching files with pattern ${pattern}:")
                    potentialFiles.each { f -> logger.info("  ${f.name}") }
                    return potentialFiles[0] // Return first match
                }
            }
        }

        logger.warning("No MAESTRO files found matching any UniProt patterns")
        return null
    }

    /**
     * Calculate Instability Heat Score (IHS) for a protein structure using MAESTRO predictions
     *
     * @param residues List of residues in the protein
     * @param pdbId The PDB ID of the protein
     * @return Calculated IHS value
     */
    private double calculateIHS(List<Residue> residues, String pdbId) {
        logger.info("\nStarting IHS calculation for PDB: ${pdbId}")

        // Clean PDB ID (remove .pdb extension if present)
        String cleanPdbId = pdbId.replaceAll(/(?i)\.pdb$/, "")
        logger.info("Cleaned PDB ID: ${cleanPdbId}")

        File maestroDir = new File(stabilityDir)
        if (!maestroDir.exists()) {
            logger.warning("Stability directory not found: ${maestroDir.absolutePath}")
            return 0.0
        }

        File maestroFile = findMaestroFile(cleanPdbId, maestroDir)
        if (!maestroFile) {
            logger.warning("No matching MAESTRO file found for ${cleanPdbId}")
            logger.warning("Searched in: ${maestroDir.absolutePath}")
            return 0.0
        }
        logger.info("Found MAESTRO file: ${maestroFile.absolutePath}")

        List<Double> significantDDGs = []
        int totalMutations = 0
        int linesProcessed = 0

        // Parse the MAESTRO file
        try {
            logger.info("Parsing MAESTRO file...")
            maestroFile.eachLine { line ->
                linesProcessed++
                if (!line.startsWith("#") && !line.trim().isEmpty()) {
                    String[] parts = line.trim().split("\\s+")
                    if (parts.size() >= 7) { // Ensure we have all expected columns
                        String mutation = parts[3]
                        if (mutation != "wildtype") {
                            totalMutations++
                            try {
                                double ddG = Math.abs(parts[6].toDouble())
                                if (ddG > 0.5) {
                                    significantDDGs.add(ddG)
                                }
                            } catch (NumberFormatException e) {
                                logger.warning("Could not parse ddG value in line ${linesProcessed}: ${line}")
                            }
                        }
                    }
                }
            }

            logger.info("MAESTRO file parsing complete")
            logger.info("Total lines processed: ${linesProcessed}")
            logger.info("Total mutations found: ${totalMutations}")
            logger.info("Significant mutations (|ΔΔG| > 0.5): ${significantDDGs.size()}")

            // Calculate IHS according to the formula: (1/|S|) * sum(|ΔΔG|)
            if (significantDDGs.isEmpty()) {
                logger.info("No significant mutations found for IHS calculation")
                return 0.0
            }

            double sumDDG = significantDDGs.sum()
            double ihs = sumDDG / significantDDGs.size()

            logger.info("IHS calculation results:")
            logger.info("  Sum of |ΔΔG| values: ${String.format("%.2f", sumDDG)}")
            logger.info("  Number of significant mutations: ${significantDDGs.size()}")
            logger.info("  Calculated IHS: ${String.format("%.2f", ihs)}")

            return ihs

        } catch (Exception e) {
            logger.warning("Error parsing MAESTRO file ${maestroFile.name}: ${e.message}")
            e.printStackTrace()
            return 0.0
        }
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
                    writer.write("Residue_Number${delimiter}Residue_Name${delimiter}Atom_Name${delimiter}Chain${delimiter}B-factor\n")
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

    private double calculateMedian(List<Map<String, Object>> bFactorData) {
        List<Double> bFactors = bFactorData.collect { it.bFactor }
        bFactors.sort()

        int size = bFactors.size()
        if (size == 0) {
            return Double.NaN
        }

        if (size % 2 == 0) {
            return (bFactors[size/2 - 1] + bFactors[size/2]) / 2.0
        } else {
            return bFactors[size/2]
        }
    }

    private void printMedianBFactor(List<Map<String, Object>> bFactorData) {
        List<Double> bFactors = bFactorData.collect { it.bFactor }
        bFactors.sort()
        double median = calculateMedian(bFactorData)

        logger.info("\n B-factor statistics:")
        logger.info(String.format("  Total atoms: %d", bFactors.size()))
        logger.info(String.format("  Minimum B-factor: %.2f", bFactors.min()))
        logger.info(String.format("  Maximum B-factor: %.2f", bFactors.max()))
        logger.info(String.format("  Median B-factor: %.2f", median))
    }

    private void writeExtendedResultsToCSV(List<Map<String, Object>> results) {
        String outputFilename = "PDB_Features.csv"

        try {
            File outputFile = new File(outputFilename)
            outputFile.withWriter { writer ->
                // Write header
                if (calculateIHS) {
                    writer.write("PDB_ID${delimiter}Median_Bfactor${delimiter}Total_SASA${delimiter}Num_Residues${delimiter}IHS\n")
                } else {
                    writer.write("PDB_ID${delimiter}Median_Bfactor${delimiter}Total_SASA${delimiter}Num_Residues\n")
                }

                // Write data
                results.each { result ->
                    if (calculateIHS) {
                        writer.write("${result.pdbId}${delimiter}${result.medianBfactor}${delimiter}${result.totalSASA}${delimiter}${result.numResidues}${delimiter}${result.IHS}\n")
                    } else {
                        writer.write("${result.pdbId}${delimiter}${result.medianBfactor}${delimiter}${result.totalSASA}${delimiter}${result.numResidues}\n")
                    }
                }
            }
        } catch (Exception e) {
            logger.severe("Error writing results to CSV: ${e.message}")
        }
    }
}