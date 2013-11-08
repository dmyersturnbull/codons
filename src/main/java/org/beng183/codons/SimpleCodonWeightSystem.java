package org.beng183.codons;

import java.util.Map;

public class SimpleCodonWeightSystem implements CodonWeightSystem {

	private Map<String,Double> map;
	
	public SimpleCodonWeightSystem() {
		
	}
	
	@Override
	public double getCodonWeight(String codon) {
		return map.get(codon);
	}

}
