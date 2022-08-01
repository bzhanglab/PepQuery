package main.java.util;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class Cloger {

	private static Cloger instance = null;
	private long start_time;

	public Logger logger;


	private Cloger(){
		start_time = System.currentTimeMillis();
		logger = LogManager.getLogger(createProDB.class.getName());
	}

	public static Cloger getInstance() {
		if (instance == null) {
			instance = new Cloger();
		}
		return instance;
	}

	public String getRunTime(){
		long ctime = System.currentTimeMillis();
		double t = 1.0*(ctime  - start_time)/1000.0/60.0;
		String out = String.format("%.2f",t);
		out = out + " min";
		return(out);
	}
}


