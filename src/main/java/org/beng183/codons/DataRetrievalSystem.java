package org.beng183.codons;

import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.scop.ScopDomain;

/**
 * A way to get genetic sequences, protein structures, and domains given a gene identifier.
 * @author dmyersturnbull
 */
public interface DataRetrievalSystem {

	/**
	 * Returns the full gene sequence for a gene identifier.
	 */
	String getGeneticSequenceFromName(String name) throws LoadException;

	/**
	 * Returns the the 3-dimensional structure of the first PDB Id mapped from the gene identifier.
	 */
	Structure getStructureFromGeneName(String name) throws LoadException;
	
	/**
	 * Returns the list of all domains mapped from the gene identifier.
	 */
	List<ScopDomain> getDomainsFromGeneName(String name) throws LoadException;
	
}
