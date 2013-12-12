package org.beng183.codons;

/**
 * The name of a species.
 * @author dmyersturnbull
 */
public enum Species {
	
	E_COLI("E. coli"), S_CEREVISIAE("S. cerevisiae");
	
	private final String name;

	private Species(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
	}
	
}
