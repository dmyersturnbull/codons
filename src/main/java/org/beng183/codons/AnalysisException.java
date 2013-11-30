package org.beng183.codons;

/**
 * An error in an evaluation of correlation. A scientific error.
 * @author dmyersturnbull
 */
@SuppressWarnings("serial")
public class AnalysisException extends Exception {

	public AnalysisException() {
		super();
	}

	public AnalysisException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) {
		super(message, cause, enableSuppression, writableStackTrace);
	}

	public AnalysisException(String message, Throwable cause) {
		super(message, cause);
	}

	public AnalysisException(String message) {
		super(message);
	}

	public AnalysisException(Throwable cause) {
		super(cause);
	}

}
