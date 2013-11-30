package org.beng183.codons;

/**
 * Essentially, a mapping of 3-letter codons to their numerical scores.
 * A <em>high weight is fast (more common)</em>, while a <em>low weight is slow (less common)</em>.
 * @author dmyersturnbull
 */
public interface CodonWeightSystem {

	/**
	 * Returns the "speed" weight of the given 3-letter codon.
	 * A <em>high weight is fast (more common)</em>, while a <em>low weight is slow (less common)</em>.
	 */
	double getCodonWeight(String codon);
	
}
