package org.beng183.codons;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.ref.WeakReference;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava3.core.sequence.DNASequence;

/**
 * A simple {@link DataRetrievalSystem} that uses Entrez to fetch genetic sequences,
 * and the <a href="http://www.ebi.ac.uk/pdbe/docs/sifts/">SIFTS</a> to map between gene names, PDB identifiers, and domain identifiers.
 * @author dmyersturnbull
 */
public class SiftsDataRetrievalSystem implements DataRetrievalSystem {

	private static final Logger logger = LogManager.getLogger(SiftsDataRetrievalSystem.class.getName());

	private Map<String, WeakReference<String>> geneToSequence = new HashMap<>();
	private Map<String, String> chainIds = new HashMap<>();
	private Map<String, String> scopIds = new HashMap<>();
	private Map<String, String> entrezToUniProt = new HashMap<>();
	private Map<String, String> entrezToRefSeq = new HashMap<>();

	private AtomCache cache;

	private ScopDatabase scop;

	public SiftsDataRetrievalSystem() {
		cache = new AtomCache();
		scop = new BerkeleyScopInstallation();
		try {
			readScopMapping();
			readChainMapping();
		} catch (IOException e) {
			throw new RuntimeException("Couldn't load SIFTS mapping", e);
		}
	}
	
	private void readChainMapping() throws IOException {
		File file = new File("src/main/resources/pdb_chain_uniprot.tsv");
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String line = "";
			while ((line = br.readLine()) != null) {
				if (line.isEmpty() || line.startsWith("#") || line.startsWith("PDB")) continue;
				String[] parts = line.split("\t");
				String pdbId = parts[0];
				String chainId = parts[1];
				String uniProtId = parts[2];
				String seqresStart = parts[3];
				String seqresEnd = parts[4];
				String pdbStart = parts[5];
				String pdbEnd = parts[6];
				String uniprotStart = parts[7];
				String uniprotEnd = parts[8];
				chainIds.put(uniProtId, pdbId + "." + chainId);
//				chainIds.put(uniProtId, pdbId + "." + chainId + "_" + pdbStart + "-" + pdbEnd);
			}
		}
	}

	private void readScopMapping() throws IOException {
		File file = new File("src/main/resources/pdb_chain_scop_uniprot.tsv");
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String line = "";
			while ((line = br.readLine()) != null) {
				if (line.isEmpty() || line.startsWith("#") || line.startsWith("PDB")) continue;
				String[] parts = line.split("\t");
				String pdbId = parts[0];
				String chainId = parts[1];
				String uniProtId = parts[2];
				int sunId = Integer.parseInt(parts[3]);
				String scopId = parts[4];
				scopIds.put(uniProtId, scopId);
			}
		}
	}

	@Override
	public String getGeneticSequenceFromName(String name) throws LoadException {

		logger.info("Retrieving sequence for gene " + name + "...");

		if (geneToSequence.containsKey(name)) {
			String ans = geneToSequence.get(name).get();
			if (ans != null) return ans;
		}

		String uniProtId = getUniProtId(name);
		String refSeqId = getRefSeqId(uniProtId);
		LinkedHashMap<String, DNASequence> sequences = SimpleDataRetrievalSystem.downloadSequences(refSeqId, name);
		
		Iterator<DNASequence> iterator = sequences.values().iterator();
		if (!iterator.hasNext()) throw new LoadException("Did not find a sequence for gene with name " + name);
		String sequence = iterator.next().getSequenceAsString();
		geneToSequence.put(name, new WeakReference<String>(sequence));
		logger.info("Found sequence of length " + sequence.length() + " for gene " + name);
		return sequence;
	}

	@Override
	public Structure getStructureFromGeneName(String name) throws LoadException {
		logger.info("Finding PDB structure for gene " + name + "...");
		String uniProtId = getUniProtId(name);
		String pdbId = chainIds.get(uniProtId);
		if (pdbId == null) throw new LoadException("PDB Id for gene " + name + " (UniProtId=" + uniProtId + ") not found");
		try {
			return cache.getStructure(pdbId);
		} catch (IOException | StructureException e) {
			throw new LoadException("Couldn't load structure from PDB Id" + pdbId, e);
		}
	}

	@Override
	public List<ScopDomain> getDomainsFromGeneName(String name) throws LoadException {
		logger.info("Finding domains for gene " + name + "...");
		String uniProtId = getUniProtId(name);
		String scopId = scopIds.get(uniProtId);
		if (scopId == null) {
			return null;
		}
		ScopDomain domain = scop.getDomainByScopID(scopId);
		List<ScopDomain> domains = new ArrayList<ScopDomain>();
		domains.add(domain);
		logger.info("Found scop Id " + scopId + " for gene " + name);
		return domains;
	}

	private String getUniProtId(String entrez) throws LoadException {
		String uniProtId = entrezToUniProt.get(entrez);
		if (uniProtId == null) {
			uniProtId = SimpleDataRetrievalSystem.getUniProtIdFromEntrezId(entrez);
			entrezToUniProt.put(entrez, uniProtId);
		}
		return uniProtId;
	}

	private String getRefSeqId(String uniProtId) throws LoadException {
		String refSeqId = entrezToRefSeq.get(uniProtId);
		if (refSeqId == null) {
			refSeqId = SimpleDataRetrievalSystem.getRefSeqIdFromUniProtId(uniProtId);
			entrezToRefSeq.put(uniProtId, refSeqId);
		}
		return refSeqId;
	}

}
