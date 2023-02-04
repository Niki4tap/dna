use dna_rs::OwnedDNA;

fn main() {
	let args = std::env::args().skip(1).collect::<Vec<String>>();

	match args[0].as_str() {
		"enc" => {
			let dna = OwnedDNA::try_from_nucleotide_str(&args[1]).unwrap();
			println!("{}", dna.try_to_byte_str().unwrap());
		},
		"dec" => {
			let dna = OwnedDNA::try_from_byte_str(&args[1]).unwrap();
			println!("{}", dna.to_nucleotide_str());
		},
		_ => {
			panic!("unrecognized command");
		}
	}
}
