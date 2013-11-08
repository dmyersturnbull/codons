package org.beng183.codons;

import java.util.List;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.scop.ScopDomain;

public interface DataRetrievalSystem {

	String getGeneticSequenceFromName(String name) throws LoadException;
	
	Structure getStructureFromGeneName(String name) throws LoadException;
	
	List<ScopDomain> getDomainsFromGeneName(String name) throws LoadException;
	
}
