package org.beng183.codons;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomPositionMap;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.ResidueRange;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.secstruc.SecStruc;
import org.biojava.bio.structure.secstruc.SecStrucGroup;
import org.biojava.bio.structure.secstruc.SecStrucState;
import org.biojava.bio.structure.secstruc.SecStrucType;

/**
 * A class to determine the correlation between codon usage bias and something else.
 * So far, this includes:
 * <ul>
 * <li>domain boundaries</li>
 * <li>α versus β SSEs</li>
 * <li>sequence length</li>
 * <li>root mean squared deviation of the residues</li>
 * </ul>
 * Methods that evaluate the correlation between codon usage bias and a binary variable output an {@link Evaluation} object to describe the correlation.
 * @author dmyersturnbull
 *
 */
public class SimpleCorrelations {

	public static void main(String[] args) throws AnalysisException, IOException {
		
		if (args.length != 1) {
			System.err.println("Usage: " + SimpleCorrelations.class.getSimpleName() + " list-of-ids-file");
			return;
		}
		
		logger.info("Initializing data retrieval system");
		DataRetrievalSystem retrieval = new SiftsDataRetrievalSystem();
		
		logger.info("Initializing codon weight system");
		CodonWeightSystem weighter = SimpleCodonWeightSystem.createForSpecies(Species.S_CEREVISIAE);
		
		logger.info("Initializing correlations-finder");
		SimpleCorrelations corr = new SimpleCorrelations(weighter, retrieval);
		
		logger.info("Parsing gene identifiers");
		corr.setNames(new File(args[0]));

		printHeader("domains");
		Evaluation domains = corr.evaluateNearDomainBoundaries(5);
		System.out.println(domains);

		printHeader("strands");
		Evaluation strands = corr.evaluateWithinBetaSheets();
		System.out.println(strands);

		printHeader("length");
		RealMatrix length = corr.correlateLength();
		System.out.println(printMatrix(length));

		printHeader("sparseness");
		RealMatrix sparseness = corr.correlateSparseness();
		System.out.println(printMatrix(sparseness));
	}

	/**
	 * An object that calculates a numerical value from a structure.
	 * This value can then be correlated against the total codon weight of the sequence.
	 * @author dmyersturnbull
	 */
	public interface StructureCalculator {
		/**
		 * Calculates the value.
		 * @param structure The structure of the whole protein
		 * @param ca An array of C-α atoms in {@code structure}
		 * @param geneticSequence The genetic sequence; should probably be used only in an auxiliary way
		 */
		double calculate(Structure structure, Atom[] ca, String geneticSequence) throws Exception;
	}

	/**
	 * A description of the correlation between codon usage bias and a binary variable.
	 * @author dmyersturnbull
	 */
	public static class Evaluation {

		// each pair is indexed by (positive, negative) (which is (true, false))
		List<Pair<DescriptiveStatistics, DescriptiveStatistics>> stats = new ArrayList<>();

		void add(DescriptiveStatistics positive, DescriptiveStatistics negative) {
			stats.add(new Pair<DescriptiveStatistics, DescriptiveStatistics>(positive, negative));
		}

		/**
		 * Returns statistics about the <em>means</em> of the <em>negative</em> proteins
		 * For example, if "positive" is β-sheet and "negative" is α-helix of n proteins,
		 * this will return statistics of the n-dimensional vector containing the
		 * <em>mean codon usage bias</em> for each protein. The mean for a protein x is the
		 * mean codon usage bias among α-helical residues.
		 */
		public DescriptiveStatistics getNegativeMeans() {
			DescriptiveStatistics means = new DescriptiveStatistics();
			for (Pair<DescriptiveStatistics, DescriptiveStatistics> pair : stats) {
				means.addValue(pair.getSecond().getMean());
			}
			return means;
		}

		/**
		 * Returns statistics about the <em>standard deviations</em> of the <em>negative</em> proteins
		 * For example, if "positive" is β-sheet and "negative" is α-helix of n proteins,
		 * this will return statistics of the n-dimensional vector containing the
		 * <em>standard deviation of codon usage bias</em> for each protein. The standard deviation
		 * for a protein x is the standard deviation of codon usage bias among α-helical residues.
		 */
		public DescriptiveStatistics getNegativeStandardDeviations() {
			DescriptiveStatistics means = new DescriptiveStatistics();
			for (Pair<DescriptiveStatistics, DescriptiveStatistics> pair : stats) {
				means.addValue(pair.getSecond().getStandardDeviation());
			}
			return means;
		}

		/**
		 * @see #getNegativeMeans()
		 */
		public DescriptiveStatistics getPositiveMeans() {
			DescriptiveStatistics means = new DescriptiveStatistics();
			for (Pair<DescriptiveStatistics, DescriptiveStatistics> pair : stats) {
				means.addValue(pair.getFirst().getMean());
			}
			return means;
		}

		/**
		 * @see #getNegativeStandardDeviations()
		 */
		public DescriptiveStatistics getPositiveStandardDeviations() {
			DescriptiveStatistics means = new DescriptiveStatistics();
			for (Pair<DescriptiveStatistics, DescriptiveStatistics> pair : stats) {
				means.addValue(pair.getFirst().getStandardDeviation());
			}
			return means;
		}
		@Override
		public String toString() {
			StringBuilder sb = new StringBuilder();
			sb.append("positive: " + getPositiveMeans().getMean() + "; ");
			sb.append("negative: " + getNegativeMeans().getMean());
			return sb.toString();
		}
	}

	private static final Logger logger = LogManager.getLogger(SimpleCorrelations.class.getName());

	private List<String> names;
	private DataRetrievalSystem retrieval;
	private CodonWeightSystem weighter;

	public SimpleCorrelations(CodonWeightSystem codonWeighter, DataRetrievalSystem retrievalSystem) {
		super();
		this.weighter = codonWeighter;
		this.retrieval = retrievalSystem;
	}

	/**
	 * Correlates codon usage bias with sequence length.
	 * @see #correlate(StructureCalculator)
	 */
	public RealMatrix correlateLength() throws AnalysisException {
		return correlate(new StructureCalculator() {
			@Override
			public double calculate(Structure structure, Atom[] ca, String geneticSequence) throws Exception {
				return ca.length;
			}
		});
	}

	/**
	 * Correlates codon usage bias with the average distance of an atom to its closest neighbor.
	 * <pre>
	 * E[ min{|i-j| : j≠i in S } : i in S]
	 * </pre>
	 * This should give an estimate of inverse density, without normalization for length.
	 * @see #correlate(StructureCalculator)
	 */
	public RealMatrix correlateSparseness() throws AnalysisException {
		return correlate(new StructureCalculator() {
			@Override
			public double calculate(Structure structure, Atom[] ca, String geneticSequence) throws Exception {
				return superpositionDistance(ca, ca);
			}
		});
	}

	/**
	 * Correlates the total codon usage bias of each structure against some other value, as calculated by {@code calculator}.
	 * @return A PearsonsCorrelation object with codon usage bias in the first column and the calculated value in the second
	 */
	public RealMatrix correlate(StructureCalculator calculator) throws AnalysisException {

		double[] calculatedValue = new double[names.size()];
		double[] weights = new double[names.size()];

		for (int x = 0; x < names.size(); x++) {

			String name = names.get(x);

			try {

				String geneticSequence = retrieval.getGeneticSequenceFromName(name);
				Structure structure = retrieval.getStructureFromGeneName(name);
				Atom[] ca = StructureTools.getAtomCAArray(structure);

				if (ca.length != geneticSequence.length() / 3) {
					logger.warn(name + " has " + ca.length + " C-α atoms but " + (geneticSequence.length()/3) + " codons");
				}

				calculatedValue[x] = calculator.calculate(structure, ca, geneticSequence);
				weights[x] = sumWeight(geneticSequence);

			} catch (Exception e) {
				logger.error(e);
				e.printStackTrace();
			}
		}

		RealMatrix matrix = new BlockRealMatrix(names.size(), 2);
		matrix.setColumn(0, calculatedValue);
		matrix.setColumn(1, weights);
		return matrix;
	}

	/**
	 * Evaluates codon usage bias between two classes:
	 * <ul>
	 * <li>Residues within a sequence distance of nResiduesAround from any domain boundary, inclusive (positive), and</li>
	 * <li>Residues outside of a sequence distance of nResiduesAround from any domain boundary (negative)</li>
	 * </ul>
	 * @param nResiduesAround
	 * @return An Evaluation object, with "positive" being near domain boundaries, and "negative" being outside
	 * @throws AnalysisException
	 */
	public Evaluation evaluateNearDomainBoundaries(int nResiduesAround) throws AnalysisException {

		if (names == null) {
			throw new IllegalArgumentException("Must set names first");
		}

		Evaluation eval = new Evaluation();

		for (String name : names) {

			DescriptiveStatistics positive = new DescriptiveStatistics();
			DescriptiveStatistics negative = new DescriptiveStatistics();

			try {

				String geneticSequence = retrieval.getGeneticSequenceFromName(name);
				Structure structure = retrieval.getStructureFromGeneName(name);
				Atom[] ca = StructureTools.getAtomCAArray(structure);

				if (ca.length != geneticSequence.length() / 3) {
					throw new AnalysisException("Failed on gene " + name + " because it has " + ca.length + " C-α atoms but " + (geneticSequence.length()/3) + " codons");
				}

				List<ScopDomain> domains = retrieval.getDomainsFromGeneName(name);
				if (domains.size() < 2) {
					logger.info("Skipping " + name + " because only " + domains.size() + " domain(s) were found");
					continue;
				}

				for (ScopDomain domain : domains) {

					AtomPositionMap atomPositionMap = new AtomPositionMap(ca);
					List<String> rangeStrings = domain.getRanges();
					List<ResidueRange> ranges = ResidueRange.parseMultiple(rangeStrings);
					ResidueNumber lastResidue = ranges.get(ranges.size() - 1).getEnd();

					for (int i = 0; i < ca.length; i++) {
						Atom atom = ca[i];
						Group group = atom.getGroup();
						double codonWeight = weighter.getCodonWeight(getCodon(geneticSequence, i));
						if (group != null) {
							ResidueNumber residueNumber = group.getResidueNumber();
							int length = atomPositionMap.calcLength(residueNumber, lastResidue);
							if (length <= nResiduesAround) {
								positive.addValue(codonWeight);
							} else {
								negative.addValue(codonWeight);
							}
						}
					}

				}

				logger.info("For gene " + name + ": bias near boundaries is " + String.format("%1$.4f", positive.getMean()) + ", bias outside is " + String.format("%1$.4f", negative.getMean()));
				eval.add(positive, negative);

			} catch (Exception e) {
				logger.error(e);
				e.printStackTrace();
			}
		}
		return eval;
	}

	/**
	 * Sets the list of gene identifiers to run on.
	 * <strong>Must be called before any evaluation.</strong>
	 * @param file A file containing a line-by-line list of gene identifiers
	 * @throws IOException 
	 * @see #setNames(List)
	 */
	public void setNames(File file) throws IOException {
		names = new ArrayList<>();
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			String line = "";
			while ((line = br.readLine()) != null) {
				if (!line.isEmpty() && !line.startsWith(";")) names.add(line);
			}
		}
	}

	/**
	 * Sets the list of gene identifiers to run on.
	 * <strong>Must be called before any evaluation.</strong>
	 * @see #setNames(File)
	 */
	public void setNames(List<String> names) {
		this.names = names;
	}

	/**
	 * Evaluates codon usage bias between two classes:
	 * <ul>
	 * <li>Residues within β-sheets (positive), and</li>
	 * <li>Residues within other residues (negative)</li>
	 * </ul>
	 * @param nResiduesAround
	 * @return An Evaluation object, with "positive" within β-sheets, and "negative" being not within β-sheets
	 * @throws AnalysisException
	 */
	public Evaluation evaluateWithinBetaSheets() throws AnalysisException {

		if (names == null) {
			throw new IllegalArgumentException("Must set names first");
		}

		Evaluation eval = new Evaluation();

		for (String name : names) {

			DescriptiveStatistics positive = new DescriptiveStatistics();
			DescriptiveStatistics negative = new DescriptiveStatistics();

			try {

				String geneticSequence = retrieval.getGeneticSequenceFromName(name);
				Structure structure = retrieval.getStructureFromGeneName(name);
				Atom[] ca = StructureTools.getAtomCAArray(structure);

				if (ca.length != geneticSequence.length() / 3) {
					throw new AnalysisException("Failed on gene " + name + " (" + structure.getPDBHeader().getIdCode() + ") because it has " + ca.length + " C-α atoms but " + (geneticSequence.length()/3) + " codons");
				}

				SecStruc ss = new SecStruc();
				try {
					ss.assign(structure);
				} catch (StructureException e) {
					throw new AnalysisException("Failed assigning secondary structure", e);
				}
				Map<ResidueNumber, SecStrucState> map = new HashMap<>();
				SecStrucGroup[] ssgs = ss.getGroups();
				for (SecStrucGroup ssg : ssgs) {
					SecStrucState state = (SecStrucState) ssg.getProperty("secstruc");
					map.put(ssg.getResidueNumber(), state);
				}

				for (int i = 0; i < ca.length; i++) {
					Atom atom = ca[i];
					Group group = atom.getGroup();
					String codon = getCodon(geneticSequence, i);
					if (codon == null) {
						logger.warn("Codon #" + i + " does not exist");
						continue;
					}
					double codonWeight = weighter.getCodonWeight(codon);
					if (group != null) {
						ResidueNumber residueNumber = group.getResidueNumber();
						SecStrucType type = map.get(residueNumber).getSecStruc();
						if (type.equals(SecStrucType.bridge) || type.equals(SecStrucType.extended)) {
							positive.addValue(codonWeight);
						} else {
							negative.addValue(codonWeight);
						}
					}
				}

				eval.add(positive, negative);

			} catch (Exception e) {
				logger.error(e);
				e.printStackTrace();
			}
		}
		return eval;
	}

	private double superpositionDistance(Atom[] ca1, Atom[] ca2) throws StructureException {
		double total = 0;
		for (int i = 0; i < ca1.length; i++) {
			for (int j = 0; j < ca2.length; j++) {
				total += Calc.getDistanceFast(ca1[i], ca2[j]);
			}
		}
		return Math.sqrt(total) / ca1.length;
	}

	/**
	 * Calculates the sum of the codon weights in the sequence.
	 */
	private double sumWeight(String geneticSequence) {
		double weightSum = 0;
		for (int i = 0; i < geneticSequence.length() / 3; i++) {
			double codonWeight = weighter.getCodonWeight(getCodon(geneticSequence, i));
			weightSum += codonWeight;
		}
		return weightSum;
	}

	/**
	 * Returns the {@code nRead}th codon’s 3-letter code.
	 */
	private String getCodon(String geneticSequence, int nReads) {
		if (geneticSequence.length() < (nReads + 1) * 3) {
			return null;
		}
		return geneticSequence.substring(nReads * 3, (nReads + 1) * 3);
	}

	private static String printMatrix(RealMatrix matrix) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < matrix.getRowDimension(); i++) {
			for (int j = 0; j < matrix.getColumnDimension(); j++) {
				sb.append((String.format("%1$6.4f", matrix.getEntry(i, j))));
			}
			sb.append(System.getProperty("line.separator"));
		}
		return sb.toString();
	}
	
	private static void printHeader(String name) {
		System.out.println(repeat(System.getProperty("line.separator"), 4));
		System.out.println(repeat("-", 80));
		System.out.println(repeat("-", 40 - name.length()/2) + name + repeat("-", 40 - name.length()/2));
		System.out.println(repeat("-", 80));
	}
	
	private static String repeat(String s, int x) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < x; i++) sb.append(s);
		return sb.toString();
	}

}
