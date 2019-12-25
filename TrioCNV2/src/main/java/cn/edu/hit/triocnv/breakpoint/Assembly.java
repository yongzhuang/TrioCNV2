package cn.edu.hit.triocnv.breakpoint;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 *
 * @author Yongzhuang Liu
 */

public class Assembly {

	private String fastqFile1;
	private String fastqFile2;
	private String outputDir;
	private int thread;

	public Assembly(String fastqFile1, String fastqFile2, String outputDir, int thread) {
		this.fastqFile1 = fastqFile1;
		this.fastqFile2 = fastqFile2;
		this.outputDir = outputDir;
		this.thread = thread;
	}

	public void runCommand() throws IOException, InterruptedException {
		String cmd1 = "sga preprocess --pe-mode 1 -o reads.pp.fastq " + this.fastqFile1 + " " + this.fastqFile2;
		String cmd2 = "sga index -a ropebwt -t " + thread + " --no-reverse reads.pp.fastq";
		String cmd3 = "sga correct -k 19 --learn -t " + thread + " -o reads.ec.fastq reads.pp.fastq";
		String cmd4 = "sga index -a ropebwt -t " + thread + " reads.ec.fastq";
		String cmd5 = "sga overlap -m 45  -t " + thread + " reads.ec.fastq";
		String cmd6 = "sga assemble -m 45 --min-branch-length 150 -o primary reads.ec.asqg.gz";
		String[] cmdarray = new String[] { cmd1, cmd2, cmd3, cmd4, cmd5, cmd6 };
		for (String command : cmdarray) {
			execute(command);
		}
	}

	private void execute(String command) throws InterruptedException, IOException {

		try {
			ProcessBuilder builder = new ProcessBuilder("/bin/bash", "-c", command);
			builder.redirectErrorStream(true);
			builder.directory(new File(outputDir));
			Process process = builder.start();
			new Thread() {
				@Override
				public void run() {
					BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
					String line = null;
					try {
						while ((line = in.readLine()) != null) {
							System.out.println(line);
						}
					} catch (IOException e) {
						e.printStackTrace();
					} finally {
						try {
							in.close();
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				}
			}.start();

			new Thread() {
				@Override
				public void run() {
					BufferedReader err = new BufferedReader(new InputStreamReader(process.getErrorStream()));
					String line = null;

					try {
						while ((line = err.readLine()) != null) {
							System.out.println(line);
						}
					} catch (IOException e) {
						e.printStackTrace();
					} finally {
						try {
							err.close();
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				}
			}.start();
			int value=process.waitFor();
			System.out.println(command+"="+value);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
