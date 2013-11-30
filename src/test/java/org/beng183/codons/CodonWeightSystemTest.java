package org.beng183.codons;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/**
 * A test for {@link CodonWeightSystem}.
 * @author dmyersturnbull
 */
@RunWith(Parameterized.class)
public class CodonWeightSystemTest {

	private static final double PRECISION = Math.pow(2, -16);

	private final CodonWeightSystem weighter;
	private final Map<String,Double> valuesToTest;

	public CodonWeightSystemTest(CodonWeightSystem weighter, Map<String,Double> valuesToTest) {
		this.weighter = weighter;
		this.valuesToTest = valuesToTest;
	}

	@Test
	public final void test() {
		for (Map.Entry<String,Double> entry : valuesToTest.entrySet()) {
			assertEquals("Wrong entry for " + entry.getKey(), entry.getValue(), weighter.getCodonWeight(entry.getKey()), PRECISION);
		}
	}

	@Parameters
	public static Collection<Object[]> getInstances() {
		List<Object[]> list = new ArrayList<>();
		Map<String, Double> valuesToTest = new LinkedHashMap<>();
		double lMean = (0.11 + 0.11 + 0.10 + 0.10 + 0.03 + 0.55) / 6; // Leucine
		valuesToTest.put("UUG", 0.11 / lMean);
		valuesToTest.put("CUU", 0.10 / lMean);
		valuesToTest.put("CUC", 0.10 / lMean);
		valuesToTest.put("UUG", 0.11 / lMean);
		valuesToTest.put("CUA", 0.03 / lMean);
		valuesToTest.put("CUG", 0.55 / lMean);
		double stopMean = (0.62 + 0.09 + 0.30) / 3;
		valuesToTest.put("UGA", 0.30 / stopMean);
		valuesToTest.put("AUG", 1.0); // Met
		list.add(new Object[] {SimpleCodonWeightSystem.createForSpecies("E. coli"), valuesToTest});
		return list;
	}
}
