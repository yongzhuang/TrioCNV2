package cn.edu.hit.triocnv.breakpoint;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import org.apache.commons.io.IOUtils;

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

		ProcessBuilder builder = new ProcessBuilder("/bin/bash", "-c", command);
		builder.redirectErrorStream(true);
		builder.directory(new File(workDir));
		Process process = builder.start();
//		try {
//			InputStream is1 = process.getInputStream();
//			InputStream is2 = process.getErrorStream();
//			new Thread() {
//				@Override
//				public void run() {
//					BufferedReader in = new BufferedReader(new InputStreamReader(is1));
//					String line = null;
//					try {
//						while ((line = in.readLine()) != null) {
//							System.out.println(line);
//						}
//					} catch (IOException e) {
//						e.printStackTrace();
//					} finally {						
//						try {
//							in.close();
//						} catch (IOException e) {
//							e.printStackTrace();
//						}
//					}
//				}
//			}.start();
//
//			new Thread() {
//				@Override
//				public void run() {
//					BufferedReader err = new BufferedReader(new InputStreamReader(is2));
//					String line = null;
//
//					try {
//						while ((line = err.readLine()) != null) {
//							System.out.println(line);
//						}
//					} catch (IOException e) {
//						e.printStackTrace();
//					} finally {
//						try {
//							err.close();
//						} catch (IOException e) {
//							e.printStackTrace();
//						}
//					}
//				}
//			}.start();
		int value = process.waitFor();
		// System.out.println(command + "=" + value);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//		finally {
		IOUtils.closeQuietly(process.getInputStream());
		IOUtils.closeQuietly(process.getOutputStream());
		IOUtils.closeQuietly(process.getErrorStream());
//		}
	}
}
