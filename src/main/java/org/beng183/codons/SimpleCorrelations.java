package org.beng183.codons;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
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
		DataRetrievalSystem retrieval = new SimpleDataRetrievalSystem();
		logger.info("Initializing codon weight system");
		CodonWeightSystem weighter = SimpleCodonWeightSystem.createForSpecies(Species.S_CEREVISIAE);
		logger.info("Initializing correlations-finder");
		SimpleCorrelations corr = new SimpleCorrelations(weighter, retrieval);
		logger.info("Parsing gene identifiers");
		corr.setNames(new File(args[0]));
		logger.info("Evaluating domain boundaries");
		Evaluation eval = corr.evaluateNearDomainBoundaries(5);
		System.out.println(eval.getPositiveMeans() + " / " + eval.getNegativeMeans());
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
	public PearsonsCorrelation correlateLength() throws AnalysisException {
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
	public PearsonsCorrelation correlateSparseness() throws AnalysisException {
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
	public PearsonsCorrelation correlate(StructureCalculator calculator) throws AnalysisException {

		double[] calculatedValue = new double[names.size()];
		double[] weights = new double[names.size()];

		for (int x = 0; x < names.size(); x++) {
			
			String name = names.get(x);
			
			try {

				String geneticSequence = retrieval.getGeneticSequenceFromName(name);
				Structure structure = retrieval.getStructureFromGeneName(name);
				Atom[] ca = StructureTools.getAtomCAArray(structure);
				calculatedValue[x] = calculator.calculate(structure, ca, geneticSequence);
				weights[x] = sumWeight(ca, geneticSequence);

			} catch (Exception e) {
				logger.error(e);
			}
		}

		RealMatrix matrix = new BlockRealMatrix(names.size(), names.size());
		matrix.setColumn(0, calculatedValue);
		matrix.setColumn(1, weights);
		PearsonsCorrelation corr = new PearsonsCorrelation(matrix);
		return corr;
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
					double codonWeight = weighter.getCodonWeight(getCodon(geneticSequence, i));
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
			}
		}
		return eval;
	}

	/**
	 * Provides a rough alignment-free metric for the similarity between two
	 * superimposed structures.
	 * 
	 * The average distance from each atom in {@code ca1} to the closest atom in {@code ca2}.
	 * 
	 * Copied from RotationOrderDetector of rcsb-symmetry.
	 *
	 * @param ca1 first structure
	 * @param ca2 second structure
	 * @return the average distance to the closest atom
	 * @throws StructureException if an error occurs finding distances between atoms
	 * @author sbliven
	 */
	private double superpositionDistance(Atom[] ca1, Atom[] ca2) throws StructureException {

		// Store the closest distance yet found
		double[] bestDist1 = new double[ca1.length];
		double[] bestDist2 = new double[ca2.length];
		Arrays.fill(bestDist1, Double.POSITIVE_INFINITY);
		Arrays.fill(bestDist2, Double.POSITIVE_INFINITY);

		for(int i=0;i<ca1.length;i++) {
			for(int j=0;j<ca2.length;j++) {
				double dist = Calc.getDistanceFast(ca1[i], ca2[j]);
				if( dist < bestDist1[i]) {
					bestDist1[i] = dist;
				}
				if( dist < bestDist2[j]) {
					bestDist2[j] = dist;
				}
			}
		}

		double total = 0;
		for(int i=0;i<ca1.length;i++) {
			total += Math.sqrt(bestDist1[i]);
		}
		for(int j=0;j<ca2.length;j++) {
			total += Math.sqrt(bestDist2[j]);
		}

		double dist = total/(ca1.length+ca2.length);
		return dist;
	}
	
	/**
	 * Calculates the sum of the codon weights in the sequence.
	 */
	private double sumWeight(Atom[] ca, String geneticSequence) {
		double weightSum = 0;
		for (int i = 0; i < ca.length; i++) {
			Atom atom = ca[i];
			Group group = atom.getGroup();
			if (group != null) {
				double codonWeight = weighter.getCodonWeight(getCodon(geneticSequence, i));
				weightSum += codonWeight;
			}
		}
		return weightSum;
	}

	/**
	 * Returns the {@code nRead}th codon’s 3-letter code.
	 */
	private String getCodon(String geneticSequence, int nReads) {
		return geneticSequence.substring(nReads * 3, (nReads + 1) * 3);
	}

}
