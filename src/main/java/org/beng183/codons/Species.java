package org.beng183.codons;

public enum Species {
	
	E_COLI("E. coli");
	
	private final String name;

	private Species(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
	}
	
}
