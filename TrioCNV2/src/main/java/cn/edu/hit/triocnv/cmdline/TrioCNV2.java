package cn.edu.hit.triocnv.cmdline;

import cn.edu.hit.triocnv.breakpoint.Alignment;
import cn.edu.hit.triocnv.breakpoint.Assembly;
import cn.edu.hit.triocnv.breakpoint.BreakPointCalling;
import cn.edu.hit.triocnv.breakpoint.ExtractSequence;
import cn.edu.hit.triocnv.breakpoint.MultiThreadBreakPointCalling;
import cn.edu.hit.triocnv.discordantreadpair.EstimateInsertSizes;
import cn.edu.hit.triocnv.discordantreadpair.TrioDRPCalling;
import cn.edu.hit.triocnv.discordantreadpair.TrioExtractDRPs;
import cn.edu.hit.triocnv.readdepth.Preprocessing;
import cn.edu.hit.triocnv.readdepth.TrioCalling;
import cn.edu.hit.triocnv.util.ParameterReader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;

/**
 *
 * @author Yongzhuang Liu
 */
public class TrioCNV2 {

	private static Logger logger = Logger.getLogger(TrioCNV2.class);

	public static void main(String[] args) throws IOException, InterruptedException {
		String usage = "\nTrioCNV2-0.1.1\n";
		usage = usage + "\nUsage: java -jar TrioCNV.jar <COMMAND> [OPTIONS]\n\n";
		usage = usage + "COMMANDS:\n" + "\tpreprocess\textract information\n"
				+ "\tcall\t\tcall CNVs in a parent-offspring trio\n" + "\trefine\t\trefine breakpoints\n";
		if (args.length > 0) {
			if (args[0].equals("preprocess") || args[0].equals("call") || args[0].equals("refine")) {
				run(args[0], args);
			} else {
				logger.error("Command is not recognized!\n" + usage);
			}
		} else {
			System.out.println(usage);
		}
	}

	private static void run(String cmd, String[] args) throws IOException, InterruptedException {
		long start = System.currentTimeMillis();
		CommandLineParser parser = new PosixParser();
		CommandLine commandLine = null;
		Options options = createOptions(cmd);
		try {
			if (options != null) {
				commandLine = parser.parse(options, args);
			}
		} catch (ParseException parseException) {
			logger.error("Invalid command line parameters!");
		}
		if (cmd.equals("preprocess")) {
			if (isValidated(commandLine, "preprocess")) {
				String reference = commandLine.getOptionValue("reference");
				String bams = commandLine.getOptionValue("bams");
				String ped = commandLine.getOptionValue("pedigree");
				String mappability = commandLine.getOptionValue("mappability");
				String output = commandLine.getOptionValue("output");
				int window = Integer.parseInt(commandLine.getOptionValue("window", "200"));
				int deviation = Integer.parseInt(commandLine.getOptionValue("deviation", "6"));
				int min_mapping_quality = Integer.parseInt(commandLine.getOptionValue("min_mapping_quality", "0"));
				File outputDir = new File(output);
				if (outputDir.exists()) {
					if (!outputDir.isDirectory()) {
						return;
					}
				} else {
					outputDir.mkdir();
				}
				File drpDir = new File(outputDir, "discordantreadpairs");
				if (drpDir.exists()) {
					if (!drpDir.isDirectory()) {
						return;
					}
				} else {
					drpDir.mkdir();
				}
				File rdDir = new File(outputDir, "readdepth");
				if (rdDir.exists()) {
					if (!rdDir.isDirectory()) {
						return;
					}
				} else {
					rdDir.mkdir();
				}
				String rdFile = rdDir.getAbsolutePath() + "/" + "ReadDepth.txt";
				logger.info("Extracting the read depth ......");
				(new Preprocessing(reference, bams, ped, mappability, rdFile, window, min_mapping_quality)).process();
				List<String> bamFileList = new ArrayList();
				BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(bams)));
				String line = null;
				while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
					bamFileList.add(line);
				}
				bufferedReader.close();
				String parameterFile = outputDir.getAbsolutePath() + "/parameters.txt";
				EstimateInsertSizes estimateInsertSizes = new EstimateInsertSizes(reference, bamFileList,
						parameterFile);
				estimateInsertSizes.estimate();
				ParameterReader parameterReader = new ParameterReader(parameterFile);
				Map<String, int[]> parameterMap = parameterReader.getParameters();
				logger.info("Extracting discordant read pairs ......");
				TrioExtractDRPs trioExtractDRPs = new TrioExtractDRPs(reference, bamFileList, drpDir.getAbsolutePath());
				trioExtractDRPs.extract(parameterMap, deviation);
			} else {
				printHelp("preprocess");
				return;
			}
		}
		if (cmd.equals("call")) {
			if (isValidated(commandLine, "call")) {
				String inputFolder = commandLine.getOptionValue("input");
				String outputFolder = commandLine.getOptionValue("output");
				String pedFile = commandLine.getOptionValue("pedigree");
				String mappabilityFile = commandLine.getOptionValue("mappability");
				String excludeFile = commandLine.getOptionValue("exclude");
				int numOfThreads = Integer.parseInt(commandLine.getOptionValue("nt", "1"));
				double minMappability = Double.parseDouble(commandLine.getOptionValue("min_mappability", "0"));
				int minDistance = Integer.parseInt(commandLine.getOptionValue("min_distance", "10000"));
				double transitionProb = Double.parseDouble(commandLine.getOptionValue("transition_prob", "0.00001"));
				double e = Double.parseDouble(commandLine.getOptionValue("mutation_rate", "0.0001"));
				double outlier = Double.parseDouble(commandLine.getOptionValue("outlier", "0.025"));
				int gcBinSize = Integer.parseInt(commandLine.getOptionValue("gc_bin_size", "1"));
				File outputDir = new File(outputFolder);
				if (outputDir.exists()) {
					if (!outputDir.isDirectory()) {
						return;
					}
				} else {
					outputDir.mkdir();
				}
				String readDepthFile = inputFolder + "/readdepth/ReadDepth.txt";
				String readDepthOutputFile = outputDir.getAbsolutePath() + "/ReadDepth.txt";
				String drpOutputFile = outputDir.getAbsolutePath() + "/DiscordantReadPairs.txt";
				String parameterFile = inputFolder + "/parameters.txt";
				String drpFolder = inputFolder + "/discordantreadpairs";
				new TrioCalling(readDepthFile, readDepthOutputFile, pedFile, mappabilityFile, minMappability,
						minDistance, transitionProb, outlier, e, gcBinSize).runMultiThreads(numOfThreads);
				ParameterReader parameterReader = new ParameterReader(parameterFile);
				Map<String, int[]> parameterMap = parameterReader.getParameters();
     			new TrioDRPCalling(drpFolder, pedFile, parameterFile, excludeFile, drpOutputFile).run();
			} else {
				printHelp("call");
				return;
			}
		}
		if (cmd.equals("refine")) {
			if (isValidated(commandLine, "refine")) {
			String reference = commandLine.getOptionValue("reference");
			String bams = commandLine.getOptionValue("bams");
			String ped = commandLine.getOptionValue("pedigree");
			String mappability = commandLine.getOptionValue("mappability");
			String outputFolder = commandLine.getOptionValue("output");
			String inputFolder = commandLine.getOptionValue("input");
			int size = Integer.parseInt(commandLine.getOptionValue("size", "400"));
			int deviation = Integer.parseInt(commandLine.getOptionValue("devitation", "3"));
			int numOfThreads = Integer.parseInt(commandLine.getOptionValue("nt", "8"));
			File outputDir = new File(outputFolder);
			if (outputDir.exists()) {
				if (!outputDir.isDirectory()) {
					return;
				}
			} else {
				outputDir.mkdir();
			}
			String drpSVFile = inputFolder + "/DiscordantReadPairs.txt";
			String readdepthSVFile = inputFolder + "/ReadDepth.txt";
			new MultiThreadBreakPointCalling(reference, bams, ped, readdepthSVFile, drpSVFile, outputFolder, size, deviation, numOfThreads).run();
			}
			else {
				printHelp("refine");
				return;
			}
		}
		long end = System.currentTimeMillis();
		logger.info("Total running time is " + (end - start) / 1000 + " seconds");
		logger.info("Done!");
	}

	private static Options createOptions(String cmd) {
		Options options = new Options();
		if (cmd.equals("preprocess")) {
			options.addOption(OptionBuilder.withLongOpt("reference").withDescription("reference genome file (required)")
					.hasArg().withArgName("FILE").create("R"));
			options.addOption(OptionBuilder.withLongOpt("bams").withDescription("bam list file (required)").hasArg()
					.withArgName("FILE").create("B"));
			options.addOption(OptionBuilder.withLongOpt("pedigree").withDescription("pedigree file (required)").hasArg()
					.withArgName("FILE").create("P"));
			options.addOption(OptionBuilder.withLongOpt("mappability").withDescription("mappability file (required)")
					.hasArg().withArgName("FILE").create("M"));
			options.addOption(OptionBuilder.withLongOpt("output").withDescription("output folder (required)").hasArg()
					.withArgName("FILE").create("O"));
			options.addOption(OptionBuilder.withLongOpt("window").withDescription("window size (optional, default 200)")
					.hasArg().withArgName("INT").create());
			options.addOption(OptionBuilder.withLongOpt("min_mapping_quality")
					.withDescription("minumum mapping quality (optional, default 0)").hasArg().withArgName("INT")
					.create());
			options.addOption(OptionBuilder.withLongOpt("deviation").withDescription("deviation (optional, default 6)")
					.hasArg().withArgName("INT").create());
			return options;
		} else if (cmd.equals("call")) {
			options.addOption(OptionBuilder.withLongOpt("input")
					.withDescription("input folder got by the preprocess step (required)").hasArg().withArgName("FILE")
					.create("I"));
			options.addOption(OptionBuilder.withLongOpt("output").withDescription("output folder (required)").hasArg()
					.withArgName("FILE").create("O"));
			options.addOption(OptionBuilder.withLongOpt("mappability").withDescription("mappability file (required)")
					.hasArg().withArgName("FILE").create("M"));
			options.addOption(OptionBuilder.withLongOpt("pedigree").withDescription("pedigree file (required)").hasArg()
					.withArgName("FILE").create("P"));
			options.addOption(OptionBuilder.withLongOpt("exclude").withDescription("exclude regions").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("min_mappability")
					.withDescription("minumum mappability(optional, default 0)").hasArg().withArgName("FLOAT")
					.create());
			options.addOption(OptionBuilder.withLongOpt("transition_prob").withDescription(
					"probability of transition between two different copy number states(optional, default 0.00001)")
					.hasArg().withArgName("FLOAT").create());
			options.addOption(OptionBuilder.withLongOpt("mutation_rate")
					.withDescription("de novo mutation rate (optional, default 0.0001)").hasArg().withArgName("FLOAT")
					.create());
			options.addOption(OptionBuilder.withLongOpt("min_distance")
					.withDescription("minumum distance to merge two adjacent CNVs (optional, default 10K)").hasArg()
					.withArgName("INT").create());
			options.addOption(OptionBuilder.withLongOpt("outlier")
					.withDescription(" the predefined percentage of outliers (optional, default 0.025)").hasArg()
					.withArgName("FLOAT").create());
			options.addOption(OptionBuilder.withLongOpt("gc_bin_size")
					.withDescription("size of gc bin by percent (optional, default 1)").hasArg().withArgName("INT")
					.create());
			options.addOption(OptionBuilder.withLongOpt("nt").withDescription("number of threads (optional, default 1)")
					.hasArg().withArgName("INT").create());
			return options;
		} else if (cmd.equals("refine")) {
			options.addOption(OptionBuilder.withLongOpt("reference").withDescription("reference genome file (required)")
					.hasArg().withArgName("FILE").create("R"));
			options.addOption(OptionBuilder.withLongOpt("bams").withDescription("bam list file (required)").hasArg()
					.withArgName("FILE").create("B"));
			options.addOption(OptionBuilder.withLongOpt("pedigree").withDescription("pedigree file (required)").hasArg()
					.withArgName("FILE").create("P"));
			options.addOption(OptionBuilder.withLongOpt("output").withDescription("output folder (required)").hasArg()
					.withArgName("FILE").create("O"));
			options.addOption(OptionBuilder.withLongOpt("input")
					.withDescription("input folder got by the preprocess step (required)").hasArg().withArgName("FILE")
					.create("I"));
			options.addOption(OptionBuilder.withLongOpt("nt").withDescription("number of threads (optional, default 1)")
					.hasArg().withArgName("INT").create());
			options.addOption(OptionBuilder.withLongOpt("size")
					.withDescription("the size of expanded breakpoint regions (optional, default 400)").hasArg()
					.withArgName("INT").create());
			options.addOption(OptionBuilder.withLongOpt("deviation").withDescription("deviation (optional, default 3)")
					.hasArg().withArgName("INT").create());
			return options;
		} else {
			return null;
		}
	}

	private static boolean isValidated(CommandLine line, String cmd) {
		boolean tag = true;
		if (cmd.equals("preprocess")) {
			if (!line.hasOption("reference") || !(new File(line.getOptionValue("reference")).isFile())) {
				logger.error("The reference genome file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("bams") || !(new File(line.getOptionValue("bams")).isFile())) {
				logger.error("The list of bam files is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("pedigree") || !(new File(line.getOptionValue("pedigree")).isFile())) {
				logger.error("The pedigree file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("mappability") || !(new File(line.getOptionValue("mappability")).isFile())) {
				logger.error("The mappability file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("output")) {
				logger.error("The output folder is not correctly specified!");
				tag = false;
			}
		}
		if (cmd.equals("call")) {
			if (!line.hasOption("input")) {
				logger.error("The input folder got by the preprocess step is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("output")) {
				logger.error("The output folder is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("pedigree") || !(new File(line.getOptionValue("pedigree")).isFile())) {
				logger.error("The pedigree file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("mappability") || !(new File(line.getOptionValue("mappability")).isFile())) {
				logger.error("The mappability file is not correctly specified!");
				tag = false;
			}
		}

		if (cmd.equals("refine")) {
			if (!line.hasOption("reference") || !(new File(line.getOptionValue("reference")).isFile())) {
				logger.error("The reference genome file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("bams") || !(new File(line.getOptionValue("bams")).isFile())) {
				logger.error("The list of bam files is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("pedigree") || !(new File(line.getOptionValue("pedigree")).isFile())) {
				logger.error("The pedigree file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("input")) {
				logger.error("The input folder got by the preprocess step is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("output")) {
				logger.error("The output folder is not correctly specified!");
				tag = false;
			}
		}
		return tag;
	}

	private static void printHelp(String command) {
		System.out.println();
		String usage1 = "usage: java -jar TrioCNV2.jar " + command + " [OPTIONS]\n\n"
				+ "-R,--reference\t<FILE>\treference genome file (required)\n"
				+ "-B,--bams\t<FILE>\tbam list file (required)\n" + "-P,--pedigree\t<FILE>\tpedigree file (required)\n"
				+ "-M,--mappability\t<FILE>\tmappability file (required)\n"
				+ "-O,--output\t<FILE>\toutput folder(required)\n"
				+ "	  --deviation\t<INT>\tdeletion insert size cutoff, median+deviation*SD(optional, default 6)\n"
				+ "   --window\t<INT>\twindow size (optional, default 200)\n"
				+ "   --min_mapping_quality\t<INT>\tminumum mapping quality (optional,default 0)\n";
		String usage2 = "usage: java -jar TrioCNV2.jar " + command + " [OPTIONS]\n\n"
				+ "-I,--input\t<FILE>\tinput folder got by the preprocess step (required)\n"
				+ "-P,--pedigree\t<FILE>\tpedigree file (required)\n"
				+ "-M,--mappability\t<FILE>\tmappability file (required)\n"
				+ "-O,--output\t<FILE>\toutput folder (required)\n" + "   --exclude\t<FILE>\texclude regions\n"
				+ "   --min_mappability\t<FLOAT>\tminumum mappability(optional, default 0)\n"
				+ "   --mutation_rate\t<FLOAT>\tde novo mutation rate (optional, default 0.0001)\n"
				+ "   --transition_prob\t<FLOAT>\tprobability of transition between two different copy number states(optional, default 0.00001)\n"
				+ "   --min_distance\t<INT>\tminumum distance to merge two adjacent CNVs (optional, default 10K)\n"
				+ "   --outlier\t<FLOAT>\tthe predefined percentage of outliers (optional, default 0.025)\n"
				+ "   --gc_bin_size\t<INT>\tsize of gc bin by percent (optional, default 1)\n"
				+ "   --nt\t<INT>\tnumber of threads (optional, default 1)\n";
		String usage3 = "usage: java -jar TrioCNV2.jar " + command + " [OPTIONS]\n\n"
				+ "-R,--reference\t<FILE>\treference genome file (required)\n"
				+ "-B,--bams\t<FILE>\tbam list file (required)\n" + "-P,--pedigree\t<FILE>\tpedigree file (required)\n"
				+ "-I,--input\t<FILE>\tinput folder got by the preprocess step (required)\n"
				+ "-O,--output\t<FILE>\toutput folder(required)\n"
				+ "	  --deviation\t<INT>\tdeletion insert size cutoff, median+deviation*SD(optional, default 3)\n"
				+ "   --size\t<INT>\tthe size of expanded breakpoint regions (optional, default 400)\n"
				+ "   --nt\t<INT>\tnumber of threads (optional, default 8)\n";
		if (command.equals("preprocess")) {
			System.out.println(usage1);
		} else if (command.equals("call")) {
			System.out.println(usage2);
		} else {
			System.out.println(usage3);
		}
	}
}
