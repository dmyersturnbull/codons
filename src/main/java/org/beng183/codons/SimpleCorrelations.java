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

public class SimpleCorrelations {

	public static void main(String[] args) throws AnalysisException {
		DataRetrievalSystem retrieval = new SimpleDataRetrievalSystem();
		CodonWeightSystem weighter = SimpleCodonWeightSystem.createForSpecies(Species.E_COLI);
		SimpleCorrelations corr = new SimpleCorrelations(weighter, retrieval);
		Evaluation eval = corr.evaluateNearDomainBoundaries(5);
		System.out.println(eval.getPositiveMeans() + " / " + eval.getNegativeMeans());
	}
	
	public static class Evaluation {
		List<Pair<DescriptiveStatistics, DescriptiveStatistics>> stats = new ArrayList<>();

		void add(DescriptiveStatistics positive, DescriptiveStatistics negative) {
			stats.add(new Pair<DescriptiveStatistics, DescriptiveStatistics>(positive, negative));
		}

		public DescriptiveStatistics getNegativeMeans() {
			DescriptiveStatistics means = new DescriptiveStatistics();
			for (Pair<DescriptiveStatistics, DescriptiveStatistics> pair : stats) {
				means.addValue(pair.getSecond().getMean());
			}
			return means;
		}

		public DescriptiveStatistics getNegativeStandardDeviations() {
			DescriptiveStatistics means = new DescriptiveStatistics();
			for (Pair<DescriptiveStatistics, DescriptiveStatistics> pair : stats) {
				means.addValue(pair.getSecond().getStandardDeviation());
			}
			return means;
		}

		public DescriptiveStatistics getPositiveMeans() {
			DescriptiveStatistics means = new DescriptiveStatistics();
			for (Pair<DescriptiveStatistics, DescriptiveStatistics> pair : stats) {
				means.addValue(pair.getFirst().getMean());
			}
			return means;
		}

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
	
	public Evaluation evaluateNearDomainBoundaries(int nResiduesAround) throws AnalysisException {

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

	public void setNames(List<String> names) {
		this.names = names;
	}

	public Evaluation evaluateWithinBetaSheets() throws AnalysisException {

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
