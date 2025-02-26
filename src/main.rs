use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use std::error::Error;

use clap::Parser;
use bio::io::fasta::Reader;
use flate2::write::GzEncoder;
use flate2::Compression;
use crossbeam_channel::unbounded;
use num_cpus;

type BoxError = Box<dyn Error + Send + Sync>;

#[derive(Clone)]
struct FastaRecord {
    header: String,
    sequence: Vec<u8>,
}

/// Enum to hold either a plain writer or a GZIP encoder.
enum OutputWriter {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<File>),
}

impl Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputWriter::Plain(writer) => writer.write(buf),
            OutputWriter::Gzip(writer) => writer.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputWriter::Plain(writer) => writer.flush(),
            OutputWriter::Gzip(writer) => writer.flush(),
        }
    }
}

#[derive(Parser)]
#[command(author, version, about)]
struct Args {
    /// Input FASTA file (can be plain, .gz, or .zst)
    input: String,
    /// Number of splits (N)
    n: usize,
    /// Use gzip compression for output if set to 1 (compression level 1)
    #[arg(short = 'z', default_value = "0")]
    gzip: u32,
}

fn main() -> Result<(), BoxError> {
    let args = Args::parse();

    // Determine the base name from the input file (e.g. "input" from "input.fa")
    let input_path = Path::new(&args.input);
    let base_name = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output");

    // Create channels for output writers (one channel per output file)
    let mut output_senders = Vec::with_capacity(args.n);
    let mut output_receivers = Vec::with_capacity(args.n);
    for _ in 0..args.n {
        let (tx, rx) = unbounded::<String>();
        output_senders.push(tx);
        output_receivers.push(rx);
    }

    // Spawn writer threads—each thread owns one output file.
    let mut writer_handles = Vec::with_capacity(args.n);
    for i in 0..args.n {
        let filename = format!(
            "{}_{}.fa{}",
            base_name,
            i,
            if args.gzip == 1 { ".gz" } else { "" }
        );
        let file = File::create(&filename).map_err(|e| Box::new(e) as BoxError)?;
        let writer = if args.gzip == 1 {
            OutputWriter::Gzip(GzEncoder::new(file, Compression::new(1)))
        } else {
            OutputWriter::Plain(BufWriter::new(file))
        };
        let rx = output_receivers[i].clone();
        let handle = std::thread::spawn(move || -> Result<(), BoxError> {
            let mut writer = writer;
            // Write all messages sent to this output channel.
            for msg in rx.iter() {
                writer
                    .write_all(msg.as_bytes())
                    .map_err(|e| Box::new(e) as BoxError)?;
            }
            // Finalize the writer.
            match writer {
                OutputWriter::Gzip(encoder) => {
                    encoder.finish().map_err(|e| Box::new(e) as BoxError)?;
                }
                OutputWriter::Plain(mut w) => {
                    w.flush().map_err(|e| Box::new(e) as BoxError)?;
                }
            }
            Ok(())
        });
        writer_handles.push(handle);
    }

    // Create a channel to send FASTA records to worker threads.
    let (record_tx, record_rx) = unbounded::<FastaRecord>();

    // Spawn a pool of worker threads (one per available CPU core).
    let num_workers = num_cpus::get();
    let mut worker_handles = Vec::with_capacity(num_workers);
    for _ in 0..num_workers {
        let record_rx = record_rx.clone();
        let output_senders = output_senders.clone();
        let n = args.n;
        let handle = std::thread::spawn(move || {
            for record in record_rx.iter() {
                // Split the record's sequence into n parts.
                let header = record.header;
                let seq = record.sequence;
                let mut parts: Vec<Vec<u8>> = vec![Vec::new(); n];
                for (i, &base) in seq.iter().enumerate() {
                    parts[i % n].push(base);
                }
                // For each non-empty part, format a FASTA record and send it.
                for (i, part) in parts.into_iter().enumerate() {
                    if !part.is_empty() {
                        let output = format!(">{}|part{}\n{}\n", header, i, String::from_utf8_lossy(&part));
                        if let Err(e) = output_senders[i].send(output) {
                            eprintln!("Failed to send output for part {}: {:?}", i, e);
                        }
                    }
                }
            }
        });
        worker_handles.push(handle);
    }

    // Open the input file with auto-decompression.
    let file = File::open(&args.input).map_err(|e| Box::new(e) as BoxError)?;
    let ext = input_path.extension().and_then(|e| e.to_str());
    let reader: Box<dyn std::io::Read> = match ext {
        Some("gz") => Box::new(flate2::read::GzDecoder::new(file)),
        Some("zst") => Box::new(zstd::stream::read::Decoder::new(file).map_err(|e| Box::new(e) as BoxError)?),
        _ => Box::new(file),
    };
    let buf_reader = BufReader::new(reader);
    let fasta_reader = Reader::new(buf_reader);

    // Read FASTA records and send them to the worker pool.
    for result in fasta_reader.records() {
        let record = result.map_err(|e| Box::new(e) as BoxError)?;
        let header = if let Some(desc) = record.desc() {
            format!("{} {}", record.id(), desc)
        } else {
            record.id().to_string()
        };
        let sequence = record.seq().to_vec();
        let fasta_record = FastaRecord { header, sequence };
        record_tx.send(fasta_record).map_err(|e| Box::new(e) as BoxError)?;
    }
    drop(record_tx); // signal that no more records will be sent

    // Wait for all worker threads to finish.
    for handle in worker_handles {
        handle.join().expect("Worker thread panicked");
    }

    // Dropping the output senders closes the channels, letting writer threads exit.
    drop(output_senders);

    // Wait for all writer threads to finish.
    for handle in writer_handles {
        let res: Result<(), BoxError> = handle.join().expect("Writer thread panicked");
        res?;
    }

    Ok(())
}
