package org.beng183.codons;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomPositionMap;
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

	public static void main(String[] args) throws AnalysisException {
		DataRetrievalSystem retrieval = new SimpleDataRetrievalSystem();
		CodonWeightSystem weighter = SimpleCodonWeightSystem.createForSpecies(Species.E_COLI);
		SimpleCorrelations corr = new SimpleCorrelations(weighter, retrieval);
		Evaluation eval = corr.evaluateNearDomainBoundaries(5);
		System.out.println(eval.getPositiveMeans() + " / " + eval.getNegativeMeans());
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
	 * Returns the {@code nRead}th codon’s 3-letter code.
	 */
	private String getCodon(String geneticSequence, int nReads) {
		return geneticSequence.substring(nReads * 3, (nReads + 1) * 3);
	}

//	public Evaluation evaluateByLength() throws AnalysisException {
//
//		Evaluation eval = new Evaluation();
//
//		DescriptiveStatistics s = new DescriptiveStatistics();
//		
//		for (String name : names) {
//			s.addValue(v);
//			eval.add(s, s);
//		}
//
//		return eval;
//	}
	
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
				Atom[] ca = StructureTools.getAllAtomArray(structure);
				List<ScopDomain> domains = retrieval.getDomainsFromGeneName(name);
				if (domains.size() < 2) continue;

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
				Atom[] ca = StructureTools.getAllAtomArray(structure);

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

}
