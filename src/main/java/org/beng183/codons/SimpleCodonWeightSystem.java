package org.beng183.codons;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * A simple {@link CodonWeightSystem} that reads a file (by default, from src/main/resources/frequencies/) to determine codon weights.
 * @author dmyersturnbull
 */
public class SimpleCodonWeightSystem implements CodonWeightSystem {

	private static final String SPECIES_DIR = "src/main/resources/frequencies/";
	private static final String SPECIES_FILENAME_EXTENSION = ".codons";
	
	public static final String E_COLI = "E. coli";

	private static final Logger logger = LogManager.getLogger(SimpleCodonWeightSystem.class.getName());

	private Map<String,Double> map;
	
	public static SimpleCodonWeightSystem createForSpecies(String speciesName) {
		String name = speciesName.toLowerCase().replaceAll("\\.", "_").replaceAll("\\s+", "");
		Map<String,Double> map = parse(new File(SPECIES_DIR + name + SPECIES_FILENAME_EXTENSION));
		return new SimpleCodonWeightSystem(map);
	}
	
	private static Map<String, Double> parse(File file) {
		
		// read the file to create a map per amino acid
		Map<String,Map<String, Double>> map = new HashMap<>();
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String aminoAcid = "";
			String line = "";
			while ((line = br.readLine()) != null) {
				if (line.trim().isEmpty()) continue;
				line = line.trim();
				if (line.startsWith(";")) continue;
				if (line.startsWith("!")) {
					aminoAcid = line.replace("!", "").trim();
					logger.info("Amino acid " + aminoAcid);
				} else {
					if (!map.containsKey(aminoAcid)) map.put(aminoAcid, new HashMap<String,Double>());
					String[] parts = line.split("\t");
					if (parts.length != 2) throw new IllegalArgumentException("Couldn't parse line " + line);
					map.get(aminoAcid).put(parts[0], Double.parseDouble(parts[1]));
				}
			}
		} catch (IOException e) {
			throw new IllegalArgumentException("Couldn't read file " + file, e);
		}
		
		// sanity checks
		if (map.size() != 21) {
			logger.warn(map.size() + " amino were found (including stop codons, should be 21)");
		}
		int x = 0;
		for (Map.Entry<String,Map<String,Double>> entry : map.entrySet()) {
			x += entry.getValue().size();
			double pSum = 0;
			for (double p : entry.getValue().values()) pSum += p;
			if (Math.abs(pSum - 1.0) > 0.01) {
				logger.warn("Codon values for " + entry.getKey() + " sum to " + pSum);
			}
		}
		if (x < 4*4*4) {
			logger.warn("Only " + x + " codons were found");
		} else if (x > 4*4*4) {
			logger.warn("Too many (" + x + ") codons were found");
		}
		
		// normalize by amino acid
		// make the mean of codon frequencies for each amino acid equal to 1, meaning no information
		Map<String,Double> finalMap = new HashMap<>();
		for (Map.Entry<String, Map<String,Double>> entry : map.entrySet()) {
			Map<String, Double> aMap = entry.getValue();
			double aMean = 0;
			for (double v : aMap.values()) aMean += v;
			aMean /=  aMap.size();
			for (Map.Entry<String, Double> anEntry : aMap.entrySet()) {
				finalMap.put(anEntry.getKey(), anEntry.getValue() / aMean);
			}
		}
		
		return finalMap;
		
	}

	public SimpleCodonWeightSystem(Map<String, Double> map) {
		this.map = map;
	}

	@Override
	public double getCodonWeight(String codon) {
		codon = codon.replaceAll("T", "U").toUpperCase();
		if (!map.containsKey(codon)) throw new IllegalArgumentException("Codon " + codon + " does not exist");
		return map.get(codon);
	}

	public static CodonWeightSystem createForSpecies(Species species) {
		return createForSpecies(species.getName());
	}

}
