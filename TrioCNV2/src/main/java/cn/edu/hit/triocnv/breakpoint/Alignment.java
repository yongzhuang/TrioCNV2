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

public class Alignment {

	private String referenceSequenceFile;
	private String workDir;
	private String contigFile;

	public Alignment(String contigFile, String referenceSequenceFile, String workDir) {
		super();
		this.referenceSequenceFile = referenceSequenceFile;
		this.workDir = workDir;
		this.contigFile = contigFile;
	}

	public void alignContig() throws InterruptedException, IOException {
		String cmd1 = "bwa index " + referenceSequenceFile;
		String cmd2 = "bwa mem -a " + referenceSequenceFile + " " + contigFile + " > primary-contigs.sam";
		String[] cmdarray = new String[] { cmd1, cmd2 };
		for (String command : cmdarray) {
			execute(command);
		}
	}

	private void execute(String command) throws InterruptedException, IOException {

		try {
			ProcessBuilder builder = new ProcessBuilder("/bin/bash", "-c", command);
			builder.directory(new File(workDir));
			builder.redirectErrorStream(true);
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
			process.waitFor();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
