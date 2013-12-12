package org.beng183.codons;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.lang.ref.WeakReference;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
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
import org.biojava3.core.sequence.io.FastaReaderHelper;

/**
 * A simple {@link DataRetrievalSystem} that uses Entrez to fetch genetic sequences,
 * and the UniProt mapping service to map between gene names, PDB identifiers, and domain identifiers.
 * @author dmyersturnbull
 */
public class SimpleDataRetrievalSystem implements DataRetrievalSystem {

	private static class ParameterNameValue {
		private final String name;
		private final String value;

		public ParameterNameValue(String name, String value) throws LoadException {
			try {
				this.name = URLEncoder.encode(name, "UTF-8");
				this.value = URLEncoder.encode(value, "UTF-8");
			} catch (UnsupportedEncodingException e) {
				throw new LoadException(e);
			}
		}
	}

	private static final Logger logger = LogManager.getLogger(SimpleDataRetrievalSystem.class.getName());

	private static final String UNIPROT_SERVER = "http://www.uniprot.org/";
	private static final String ENTREZ_EFETCH_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=";

	private static String run(String tool, ParameterNameValue[] params) throws LoadException {

		StringBuilder message = new StringBuilder("[");
		for (int i = 0; i < params.length; i++) {
			message.append(params[i].name + "=" + params[i].value);
			if (i < params.length - 1) message.append("; ");
		}
		message.append("]");
		logger.debug("Query " + tool + " with " + message);

		StringBuilder locationBuilder = new StringBuilder(UNIPROT_SERVER + tool + "/?");
		for (int i = 0; i < params.length; i++) {
			if (i > 0) locationBuilder.append('&');
			locationBuilder.append(params[i].name).append('=').append(params[i].value);
		}

		String location = locationBuilder.toString();
		URL url;
		try {
			url = new URL(location);
		} catch (MalformedURLException e) {
			throw new LoadException(e);
		}

		HttpURLConnection conn = null;

		try {

			// submit
			conn = (HttpURLConnection) url.openConnection();
			HttpURLConnection.setFollowRedirects(true);
			conn.setDoInput(true);
			logger.trace("Connecting...");
			conn.connect();

			// wait for a response
			int status = conn.getResponseCode();
			while (true) {
				int wait = 0;
				String header = conn.getHeaderField("Retry-After");
				if (header != null) wait = Integer.valueOf(header);
				if (wait == 0) break;
				logger.trace("Waiting (" + wait + ")...");
				conn.disconnect();
				Thread.sleep(wait * 1000);
				conn = (HttpURLConnection) new URL(location).openConnection();
				conn.setDoInput(true);
				conn.connect();
				status = conn.getResponseCode();
			}

			if (status == HttpURLConnection.HTTP_OK) {
				logger.trace("Received HTTP 200.");
				InputStream stream = conn.getInputStream();
				URLConnection.guessContentTypeFromStream(stream);
				StringBuilder sb = new StringBuilder();
				String line = "";
				try (BufferedReader br = new BufferedReader(new InputStreamReader(stream))) {
					while ((line = br.readLine()) != null) {
						sb.append(line + System.getProperty("line.separator"));
					}
				}

				return sb.toString();

			} else {
				throw new LoadException("Mapping failed. Got response " + conn.getResponseMessage());
			}

		} catch (IOException | InterruptedException e) {
			throw new LoadException(e);
		} finally {
			if (conn != null) conn.disconnect();
		}
	}

	private AtomCache cache;

	private ScopDatabase scop;

	// TODO these should really be soft caches
	private Map<String, String> geneToPdbId = new HashMap<>();
	private Map<String, WeakReference<String>> geneToSequence = new HashMap<>();

	public SimpleDataRetrievalSystem() {
		cache = new AtomCache();
		scop = new BerkeleyScopInstallation();
	}

	@Override
	public List<ScopDomain> getDomainsFromGeneName(String name) throws LoadException {
		logger.info("Finding domains for gene " + name + "...");
		String pdbId = getPdbIdFromEntrezId(name);
		if (pdbId == null) throw new LoadException("Gene " + name + " not found");
		List<ScopDomain> domains = scop.getDomainsForPDB(pdbId);
		StringBuilder message = new StringBuilder();
		for (int i = 0; i < domains.size(); i++) {
			message.append(domains.get(i).getScopId());
			if (i < domains.size() - 1) message.append(", ");
		}
		logger.info("Found scop Ids " + message.toString() + " for gene " + name);
		return domains;
	}

	@Override
	public String getGeneticSequenceFromName(String name) throws LoadException {

		logger.info("Retrieving sequence for gene " + name + "...");

		if (geneToSequence.containsKey(name)) {
			String ans = geneToSequence.get(name).get();
			if (ans != null) return ans;
		}

		String uniProtIdData = run("mapping", new ParameterNameValue[] { new ParameterNameValue("from", "P_ENTREZGENEID"),
				new ParameterNameValue("to", "ID"), new ParameterNameValue("format", "tab"),
				new ParameterNameValue("query", name) });
		String uniProtId = parse(uniProtIdData);
		String refseqData = run("mapping", new ParameterNameValue[] { new ParameterNameValue("from", "ACC+ID"),
				new ParameterNameValue("to", "REFSEQ_NT_ID"), new ParameterNameValue("format", "tab"),
				new ParameterNameValue("query", uniProtId) });
		String id = parse(refseqData);

		// http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=34577062,24475906&rettype=fasta&retmode=text
		URL url;
		try {
			url = new URL(ENTREZ_EFETCH_URL + id);
		} catch (MalformedURLException e) {
			throw new LoadException("Couldn't query for gene sequence for the gene with name " + name, e);
		}

		LinkedHashMap<String, DNASequence> sequences = new LinkedHashMap<>();
		InputStream stream;
		try {
			stream = url.openStream();
		} catch (IOException e) {
			throw new LoadException("Couldn't load URL " + url + " for gene with name " + name, e);
		}

		try {
			sequences = FastaReaderHelper.readFastaDNASequence(stream);
		} catch (Exception e) {
			throw new LoadException("Couldn't load gene sequence for the gene with name " + name, e);
		}
		Iterator<DNASequence> iterator = sequences.values().iterator();
		if (!iterator.hasNext()) throw new LoadException("Did not find a sequence for gene with name " + name);
		String sequence = iterator.next().getSequenceAsString();
		geneToSequence.put(name, new WeakReference<String>(sequence));
		logger.info("Found sequence of length " + sequence.length() + " for gene " + name);
		return sequence;
	}

	private String getPdbIdFromEntrezId(String name) throws LoadException {
		if (geneToPdbId.get(name) != null) return geneToPdbId.get(name);
		String uniProtIdData = run("mapping", new ParameterNameValue[] { new ParameterNameValue("from", "P_ENTREZGENEID"),
				new ParameterNameValue("to", "ID"), new ParameterNameValue("format", "tab"),
				new ParameterNameValue("query", name) });
		String uniProtId;
		try {
			uniProtId = parse(uniProtIdData);
		} catch (LoadException e) {
			throw new LoadException("Couldn't get UniProt Id from [" + uniProtIdData + "] for gene " + name);
		}
		String pdbIdData = run("mapping", new ParameterNameValue[] { new ParameterNameValue("from", "ACC+ID"),
				new ParameterNameValue("to", "PDB_ID"), new ParameterNameValue("format", "tab"),
				new ParameterNameValue("query", uniProtId) });
		String pdbId;
		try {
			pdbId = parse(pdbIdData);
		} catch (LoadException e) {
			throw new LoadException("Couldn't get PDB Id from [" + pdbIdData + "] for gene " + name + " (from UniProtId " + uniProtId + ")");
		}
		geneToPdbId.put(name, pdbId);
		return pdbId;
	}

	private String parse(String data) throws LoadException {
		try {
			return data.split(System.getProperty("line.separator"))[1].split("\t")[1];
		} catch (RuntimeException e) {
			throw new LoadException("Couldn't parse output " + data, e);
		}
	}

	@Override
	public Structure getStructureFromGeneName(String name) throws LoadException {
		logger.info("Finding PDB structure for gene " + name + "...");
		String pdbId = getPdbIdFromEntrezId(name);
		if (pdbId == null) throw new LoadException("Gene " + name + " not found");
		Structure structure = null;
		try {
			structure = cache.getStructure(pdbId);
		} catch (IOException | StructureException e) {
			throw new LoadException("Couldn't load structure from PDB Id" + pdbId, e);
		}
		logger.info("Found PDB Id " + structure.getPDBCode() + " for gene " + name);
		return structure;
	}

}
