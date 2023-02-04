#![feature(iter_array_chunks, lint_reasons, array_try_map)]
#![warn(clippy::pedantic)]
#![allow(clippy::transmute_ptr_to_ptr, reason = "imo transmutes are more readable")]

use core::{
	borrow::Borrow,
	fmt::{self, Debug, Display, Formatter},
	mem::transmute
};
use std::{
	num::ParseIntError,
	ops::{Deref, Index},
	str::from_utf8_unchecked
};

#[macro_export]
macro_rules! dna {
	($($e:expr)*, $n:expr) => {
		OwnedDNA::new(vec![$($e),*], $n)
	};

	($($e:expr)*) => {
		{
			let v = vec![$($e),*];
			let l = v.len() * 4;
			OwnedDNA::new(v, l)
		}
	};
}

#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Nucleotide {
	Adenine  = 0,
	Cytosine = 1,
	Guanine  = 2,
	Thymine  = 3
}

impl TryFrom<char> for Nucleotide {
	type Error = ();

	fn try_from(value: char) -> Result<Self, Self::Error> {
		match value {
			'A' | 'a' => Ok(Self::A),
			'T' | 't' => Ok(Self::T),
			'C' | 'c' => Ok(Self::C),
			'G' | 'g' => Ok(Self::G),
			_ => Err(())
		}
	}
}

impl From<Nucleotide> for char {
	fn from(value: Nucleotide) -> Self {
		match value {
			Nucleotide::Adenine => 'A',
			Nucleotide::Cytosine => 'C',
			Nucleotide::Guanine => 'G',
			Nucleotide::Thymine => 'T'
		}
	}
}

impl Nucleotide {
	pub const A: Self = Self::Adenine;
	pub const C: Self = Self::Cytosine;
	pub const G: Self = Self::Guanine;
	pub const T: Self = Self::Thymine;

	#[must_use]
	pub fn as_u8(&self) -> u8 { (*self).into() }
}

impl From<Nucleotide> for u8 {
	fn from(value: Nucleotide) -> Self { unsafe { transmute(value) } }
}

impl Debug for Nucleotide {
	fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
		match self {
			Self::Adenine => write!(f, "Nucleotide::Adenine"),
			Self::Thymine => write!(f, "Nucleotide::Thymine"),
			Self::Guanine => write!(f, "Nucleotide::Guanine"),
			Self::Cytosine => write!(f, "Nucleotide::Cytosine")
		}
	}
}

impl Display for Nucleotide {
	fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
		match self {
			Self::Adenine => write!(f, "Adenine"),
			Self::Thymine => write!(f, "Thymine"),
			Self::Guanine => write!(f, "Guanine"),
			Self::Cytosine => write!(f, "Cytosine")
		}
	}
}

impl From<NucleotidePack> for Nucleotide {
	fn from(pack: NucleotidePack) -> Self {
		let pack = pack.inner;
		unsafe { transmute(pack >> 6) }
	}
}

macro_rules! array_impl {
	(from $n:expr => $i:ident: ($($e:expr),+)) => {
		impl From<NucleotidePack> for [Nucleotide; $n] {
			fn from($i: NucleotidePack) -> Self {
				let $i = $i.inner;
				unsafe { [$(transmute($e)),+] }
			}
		}
	};

	(to $n:expr => $i:ident: ($($e:expr),+)) => {
		impl From<[Nucleotide; $n]> for NucleotidePack {
			fn from($i: [Nucleotide; $n]) -> Self {
				{$($e)|+}.into()
			}
		}
	};

	($t:tt $i:ident -> $($n:expr => ($($e:expr),+)),+) => {
		$(
			array_impl!($t $n => $i: ($($e),+));
		)+
	};
}

array_impl!(
	from pack ->
		1 => (pack >> 6),
		2 => (pack >> 6, (pack >> 4 & 0b00_00_00_11)),
		3 => (pack >> 6, (pack >> 4 & 0b00_00_00_11), (pack >> 2 & 0b00_00_00_11)),
		4 => (pack >> 6, (pack >> 4 & 0b00_00_00_11), (pack >> 2 & 0b00_00_00_11), (pack & 0b00_00_00_11))
);

array_impl!(
	to pack ->
		1 => (pack[0].as_u8() << 6),
		2 => (pack[0].as_u8() << 6, pack[1].as_u8() << 4),
		3 => (pack[0].as_u8() << 6, pack[1].as_u8() << 4, pack[2].as_u8() << 2),
		4 => (pack[0].as_u8() << 6, pack[1].as_u8() << 4, pack[2].as_u8() << 2, pack[3].as_u8())
);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(transparent)]
pub struct NucleotidePack {
	pub inner: u8
}

impl NucleotidePack {
	#[must_use]
	pub fn to_array(&self) -> [Nucleotide; 4] { (*self).into() }
}

impl Display for NucleotidePack {
	fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
		write!(f, "NucleotidePack(")?;
		f.debug_list().entries(<[Nucleotide; 4]>::from(*self)).finish()?;
		write!(f, ")")
	}
}

impl From<&[Nucleotide]> for NucleotidePack {
	fn from(pack: &[Nucleotide]) -> Self {
		assert!(
			pack.len() < 5,
			"`NucleotidePack` can contain maximum of 4 nucleotides, tried to pack {}",
			pack.len()
		);

		Self::from(match pack.len() {
			0 => 0b00_00_00_00,
			1 => pack[0].as_u8() << 6,
			2 => pack[0].as_u8() << 6 | pack[1].as_u8() << 4,
			3 => pack[0].as_u8() << 6 | pack[1].as_u8() << 4 | pack[2].as_u8() << 2,
			4 => pack[0].as_u8() << 6 | pack[1].as_u8() << 4 | pack[2].as_u8() << 2 | pack[3].as_u8(),
			_ => unreachable!()
		})
	}
}

impl From<u8> for NucleotidePack {
	fn from(value: u8) -> Self { Self { inner: value } }
}

impl Index<u8> for NucleotidePack {
	type Output = Nucleotide;

	fn index(&self, index: u8) -> &Self::Output {
		assert!(
			index < 4,
			"index out of bounds: `NucleotidePack` contains only 4 values, but the index is {index}"
		);

		let out = match index {
			0 => self.inner >> 6,
			1 => (self.inner >> 4) & 0b00_00_00_11,
			2 => (self.inner >> 2) & 0b00_00_00_11,
			3 => self.inner & 0b00_00_00_11,
			_ => unreachable!()
		};

		unsafe { transmute(&out) }
	}
}

pub struct NucleotideSlice {
	bytes: [u8]
}

pub type DNA = NucleotideSlice;

impl DNA {
	#[must_use]
	pub fn new<T: AsRef<Self> + ?Sized>(inner: &T) -> &Self { inner.as_ref() }

	#[must_use]
	pub fn iter(&self, size: usize) -> NucleotideIterator<'_> {
		NucleotideIterator {
			size,
			offset: 0,
			dna: self
		}
	}

	#[must_use]
	pub fn pack_iter(&self) -> NucleotidePackIterator<'_> {
		NucleotidePackIterator {
			offset: 0,
			dna:    self
		}
	}
}

impl AsRef<DNA> for [u8] {
	fn as_ref(&self) -> &DNA { unsafe { transmute(self) } }
}

impl<'a> From<&'a [u8]> for &'a DNA {
	fn from(value: &'a [u8]) -> Self { value.as_ref() }
}

impl Index<(usize, usize)> for DNA {
	type Output = Nucleotide;

	fn index(&self, (index, size): (usize, usize)) -> &Self::Output {
		assert!(
			size > index,
			"index out of bounds: the len is {size}, but the index is {index}"
		);

		let pack_idx = index / 4;

		let pack: &NucleotidePack = unsafe { transmute(self.bytes.index(pack_idx)) };

		#[allow(clippy::cast_possible_truncation, reason = "this should never truncate")]
		&pack[(index % 4) as u8]
	}
}

pub struct NucleotideIterator<'a> {
	size:   usize,
	offset: usize,
	dna:    &'a DNA
}

impl Iterator for NucleotideIterator<'_> {
	type Item = Nucleotide;

	fn next(&mut self) -> Option<Self::Item> {
		if self.size == self.offset {
			return None
		}

		let tmp = self.dna[(self.offset, self.size)];

		self.offset += 1;

		Some(tmp)
	}
}
pub struct NucleotidePackIterator<'a> {
	offset: usize,
	dna:    &'a DNA
}

impl Iterator for NucleotidePackIterator<'_> {
	type Item = NucleotidePack;

	fn next(&mut self) -> Option<Self::Item> {
		if self.offset == self.dna.bytes.len() {
			return None
		}

		let tmp = self.dna.bytes[self.offset];

		self.offset += 1;

		Some(NucleotidePack::from(tmp))
	}
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedDNA {
	bytes: Vec<u8>,
	size:  usize
}

impl OwnedDNA {
	#[must_use]
	pub fn new(bytes: Vec<u8>, size: usize) -> Self { Self { bytes, size } }

	#[must_use]
	pub fn iter(&self) -> NucleotideIterator<'_> {
		NucleotideIterator {
			size:   self.size,
			offset: 0,
			dna:    self
		}
	}

	#[must_use]
	pub fn pack_iter(&self) -> NucleotidePackIterator<'_> {
		NucleotidePackIterator {
			offset: 0,
			dna:    self
		}
	}

	#[must_use]
	pub fn as_slice(&self) -> &NucleotideSlice { self }

	#[must_use]
	pub fn try_to_byte_str(&self) -> Result<String, ()> {
		if self.size > u8::MAX as usize {
			return Err(())
		}

		let mut s = String::new();

		s.push_str(&format!("{:X}", self.size));

		for pack in self.pack_iter() {
			s.push_str(&format!("{:X}", pack.inner))
		}

		Ok(s)
	}

	#[must_use]
	pub fn try_from_byte_str(s: &str) -> Result<Self, DNAParseError> {
		if s.len() < 2 {
			return Err(DNAParseError::Empty)
		}

		let len = &s[0..1];
		let s = &s[1..];

		Ok(OwnedDNA::new(
			s.as_bytes()
				.chunks(2)
				.map(|s| unsafe { from_utf8_unchecked(s) })
				.filter(|s| !s.is_empty())
				.map(|s| u8::from_str_radix(s, 16))
				.collect::<Result<Vec<u8>, _>>()?,
			usize::from_str_radix(len, 16)?
		))
	}

	#[must_use]
	pub fn try_from_nucleotide_str(s: &str) -> Option<Self> {
		let mut v = vec![];

		let len = s.len();

		let mut iter = s.chars().array_chunks::<4>();

		for chars in iter.by_ref() {
			v.push(NucleotidePack::from(chars.try_map(|c| Nucleotide::try_from(c).ok())?).inner);
		}

		if let Some(remainder) = iter.into_remainder() {
			let remainder = remainder.as_slice();

			if remainder.is_empty() {
				return Some(Self::new(v, len))
			}

			v.push(
				NucleotidePack::from(
					remainder
						.iter()
						.map(|c| Nucleotide::try_from(*c).ok())
						.collect::<Option<Vec<_>>>()?
						.as_slice()
				)
				.inner
			);
		}

		Some(Self::new(v, len))
	}

	#[must_use]
	pub fn to_nucleotide_str(&self) -> String {
		let mut s = String::new();

		for nucleotide in self.iter() {
			s.push(nucleotide.into());
		}

		s
	}
}

impl Borrow<DNA> for OwnedDNA {
	fn borrow(&self) -> &DNA { self.bytes.as_slice().into() }
}

impl Deref for OwnedDNA {
	type Target = DNA;

	fn deref(&self) -> &Self::Target { self.borrow() }
}

impl FromIterator<Nucleotide> for OwnedDNA {
	fn from_iter<T: IntoIterator<Item = Nucleotide>>(iter: T) -> Self {
		let mut v = vec![];
		let mut size = 0;
		let mut iter = iter.into_iter().array_chunks::<4>();

		for array in iter.by_ref() {
			size += 4;
			v.push(NucleotidePack::from(array).inner);
		}

		if let Some(remainder) = iter.into_remainder() {
			let remainder = remainder.as_slice();

			if remainder.is_empty() {
				return Self::new(v, size)
			}

			size += remainder.len();
			v.push(NucleotidePack::from(remainder).inner);
		}

		Self::new(v, size)
	}
}

#[derive(Debug, Clone)]
pub enum DNAParseError {
	Empty,
	ParseIntError(ParseIntError)
}

impl From<ParseIntError> for DNAParseError {
	fn from(value: ParseIntError) -> Self { Self::ParseIntError(value) }
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn pack() {
		assert_eq!(NucleotidePack::from(0b00_01_10_11).to_array(), [
			Nucleotide::Adenine,
			Nucleotide::Cytosine,
			Nucleotide::Guanine,
			Nucleotide::Thymine
		]);

		assert_eq!(
			NucleotidePack::from([
				Nucleotide::Adenine,
				Nucleotide::Cytosine,
				Nucleotide::Guanine,
				Nucleotide::Thymine
			]),
			NucleotidePack::from(0b00_01_10_11)
		);
	}

	#[test]
	fn pack_index() {
		let pack = NucleotidePack::from(0b_00_01_10_11);

		assert_eq!([pack[0], pack[1], pack[2], pack[3]], pack.to_array())
	}

	#[test]
	fn iter() {
		let dna = dna![0b00_01_10_11 0b_11_10_01_00];

		assert_eq!(dna.iter().collect::<OwnedDNA>(), dna);

		let dna = dna![0b00_01_10_11 0b_11_10_01_00 0b_11_00_00_00, 9];

		assert_eq!(dna.iter().collect::<OwnedDNA>(), dna);
	}

	#[test]
	fn slice_index() {
		let dna = dna![0b00_00_00_11];

		assert_eq!(dna.as_slice()[(3, 4)], Nucleotide::Thymine);
	}

	#[test]
	fn from_byte_str() {
		let input = "33";

		assert_eq!(OwnedDNA::try_from_byte_str(input).unwrap(), dna![
			0b00_00_00_11,
			3
		]);

		let input = "91BE4C0";

		assert_eq!(
			OwnedDNA::try_from_byte_str(input).unwrap(),
			dna![0b00_01_10_11 0b_11_10_01_00 0b_11_00_00_00, 9]
		);
	}

	#[test]
	fn to_byte_str() {
		assert_eq!(
			dna![0b00_01_11_10 0b10_00_00_01 0b00_10_00_00, 10]
				.try_to_byte_str()
				.unwrap(),
			"A1E8120"
		);
	}

	#[test]
	fn from_nucleotide_str() {
		assert_eq!(
			OwnedDNA::try_from_nucleotide_str("ACTGGAACAG").unwrap(),
			dna![0b00_01_11_10 0b10_00_00_01 0b00_10_00_00, 10]
		);
	}

	#[test]
	fn to_nucleotide_str() {
		assert_eq!(
			dna![0b00_01_11_10 0b10_00_00_01 0b00_10_00_00, 10].to_nucleotide_str(),
			"ACTGGAACAG"
		);
	}
}
