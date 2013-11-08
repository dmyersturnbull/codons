package org.beng183.codons;

import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;

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

	private static String run(String tool, ParameterNameValue[] params) throws LoadException {

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
			logger.debug("Connecting...");
			conn.connect();

			// wait for a response
			int status = conn.getResponseCode();
			while (true) {
				int wait = 0;
				String header = conn.getHeaderField("Retry-After");
				if (header != null) wait = Integer.valueOf(header);
				if (wait == 0) break;
				logger.debug("Waiting (" + wait + ")...");
				conn.disconnect();
				Thread.sleep(wait * 1000);
				conn = (HttpURLConnection) new URL(location).openConnection();
				conn.setDoInput(true);
				conn.connect();
				status = conn.getResponseCode();
			}

			if (status == HttpURLConnection.HTTP_OK) {
				logger.debug("Received HTTP 200.");
				InputStream reader = conn.getInputStream();
				URLConnection.guessContentTypeFromStream(reader);
				StringBuilder builder = new StringBuilder();
				int a = 0;
				while ((a = reader.read()) != -1) {
					builder.append((char) a);
				}

				return builder.toString();

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

	public SimpleDataRetrievalSystem() {
		cache = new AtomCache();
		scop = new BerkeleyScopInstallation();
	}

	@Override
	public List<ScopDomain> getDomainsFromGeneName(String name) throws LoadException {
		String pdbId = getPdbIdFromEntrezId(name);
		if (pdbId == null) throw new LoadException("Gene " + name + " not found");
		return scop.getDomainsForPDB(pdbId);
	}

	@Override
	public String getGeneticSequenceFromName(String name) throws LoadException {
		// TODO Auto-generated method stub
		// http://www.uniprot.org/uniprot/P12345.fasta
		return null;
	}

	private String getPdbIdFromEntrezId(String name) throws LoadException {
		// TODO Should use a single query followed by a map
		// downloading this many times is extremely inefficient
		return run("mapping", new ParameterNameValue[] { new ParameterNameValue("from", "P_ENTREZGENEID"),
				new ParameterNameValue("to", "PDB_ID"), new ParameterNameValue("format", "tab"),
				new ParameterNameValue("query", name) });
	}

	@Override
	public Structure getStructureFromGeneName(String name) throws LoadException {
		String pdbId = getPdbIdFromEntrezId(name);
		if (pdbId == null) throw new LoadException("Gene " + name + " not found");
		try {
			return cache.getStructure(pdbId);
		} catch (IOException | StructureException e) {
			throw new LoadException("Couldn't load structure from PDB Id" + pdbId, e);
		}
	}

}
