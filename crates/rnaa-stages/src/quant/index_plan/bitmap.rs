//! Safe cardinality scanning for CRoaring native/portable and Bifrost local tags.
use super::reader::{Reader, add};
use anyhow::{Context, Result, bail};

#[derive(Default, Clone, Copy)]
pub(super) struct Bitmap {
    pub values: u64,
    pub maximum: Option<u32>,
}

pub(super) fn native(reader: &mut Reader<'_>, bytes: u64) -> Result<Bitmap> {
    reader.section(bytes, |r| match r.u8()? {
        1 => {
            let count = u64::from(r.u32()?);
            r.fits(count, 4)?;
            let mut maximum = None;
            for _ in 0..count {
                let value = r.u32()?;
                maximum = Some(maximum.map_or(value, |old: u32| old.max(value)));
            }
            // Native arrays may be unsorted/duplicated. Their serialized count
            // is a conservative decoded upper bound; actual loading normalizes.
            Ok(Bitmap {
                values: count,
                maximum,
            })
        }
        2 => portable(r),
        _ => bail!("invalid native index bitmap tag"),
    })
}

pub(super) fn portable(r: &mut Reader<'_>) -> Result<Bitmap> {
    let start = r.pos;
    let cookie = r.u32()?;
    let runs = cookie & 0xffff == 12347;
    let count = if runs {
        u64::from(cookie >> 16) + 1
    } else if cookie == 12346 {
        u64::from(r.u32()?)
    } else {
        bail!("invalid portable index bitmap cookie");
    };
    if count > 65536 {
        bail!("portable bitmap has too many containers");
    }
    let mut flags = [0_u8; 8192];
    if runs {
        for flag in flags.iter_mut().take(count.div_ceil(8) as usize) {
            *flag = r.u8()?;
        }
    }
    let descriptions = r.pos;
    r.skip(r.fits(count, 4)?)?;
    let offsets = r.pos;
    let has_offsets = !runs || count >= 4;
    if has_offsets {
        r.skip(r.fits(count, 4)?)?;
    }
    let mut payload = r.pos;
    let mut result = Bitmap::default();
    let mut previous_key = None;
    for i in 0..count {
        r.seek(descriptions + 4 * i)?;
        let key = u32::from(r.u16()?);
        let card = u64::from(r.u16()?) + 1;
        if previous_key.is_some_and(|old| old >= key) {
            bail!("bitmap container keys are not increasing");
        }
        previous_key = Some(key);
        if has_offsets {
            r.seek(offsets + 4 * i)?;
            if u64::from(r.u32()?) != payload - start {
                bail!("bitmap container offset differs from its bounded payload");
            }
        }
        r.seek(payload)?;
        let run = flags[(i / 8) as usize] & (1 << (i % 8)) != 0;
        let maximum = if run {
            let n = u64::from(r.u16()?);
            r.fits(n, 4)?;
            let mut actual = 0;
            let mut previous = None;
            for _ in 0..n {
                let first = u32::from(r.u16()?);
                let last = first + u32::from(r.u16()?);
                if last > 65535 || previous.is_some_and(|old| old >= first) {
                    bail!("invalid bitmap run range");
                }
                add(&mut actual, u64::from(last - first) + 1)?;
                previous = Some(last);
            }
            if actual != card {
                bail!("bitmap run cardinality differs from its header");
            }
            previous.unwrap_or(0)
        } else if card <= 4096 {
            r.fits(card, 2)?;
            let mut previous = None;
            for _ in 0..card {
                let value = u32::from(r.u16()?);
                if previous.is_some_and(|old| old >= value) {
                    bail!("bitmap array values are not increasing");
                }
                previous = Some(value);
            }
            previous.unwrap_or(0)
        } else {
            let mut actual = 0_u64;
            let mut maximum = 0;
            for word in 0..1024_u32 {
                let bits = r.u64()?;
                actual += u64::from(bits.count_ones());
                if bits != 0 {
                    maximum = word * 64 + 63 - bits.leading_zeros();
                }
            }
            if actual != card {
                bail!("bitmap bitset cardinality differs from its header");
            }
            maximum
        };
        add(&mut result.values, card)?;
        result.maximum = Some((key << 16) | maximum);
        payload = r.pos;
    }
    r.seek(payload)?;
    Ok(result)
}

pub(super) fn compatible(r: &mut Reader<'_>, bytes: u64) -> Result<Bitmap> {
    r.section(bytes, |r| {
        let position = r.pos;
        let first = r.u8()?;
        r.seek(position)?;
        if first == 1 || first == 2 {
            native(r, bytes)
        } else {
            portable(r)
        }
    })
}

pub(super) fn bitcontainer(r: &mut Reader<'_>) -> Result<Bitmap> {
    let tag = r.u64()?;
    match tag & 7 {
        1 => Ok(Bitmap {
            values: u64::from((tag >> 3).count_ones()),
            maximum: None,
        }),
        2 => Ok(Bitmap {
            values: 1,
            maximum: Some((tag >> 3) as u32),
        }),
        3 => compatible(r, tag >> 3),
        0 => {
            let header = r.u16()?;
            let words = usize::from((header & 0xfff8) >> 3);
            if words == 0 {
                return Ok(Bitmap::default());
            }
            if words < 3 {
                bail!("truncated tiny bitmap header");
            }
            let mut data = [0_u16; 8191];
            data[0] = header;
            for value in data.iter_mut().take(words).skip(1) {
                *value = r.u16()?;
            }
            let card = usize::from(data[1]);
            match header & 6 {
                0 if data[3..words]
                    .iter()
                    .map(|v| v.count_ones() as usize)
                    .sum::<usize>()
                    >= card => {}
                2 if card + 3 <= words => {}
                4 => {
                    let mut actual = 0_usize;
                    for pair in data[3..words].as_chunks::<2>().0 {
                        if pair[0] > pair[1] {
                            bail!("invalid tiny bitmap run");
                        }
                        actual += usize::from(pair[1] - pair[0]) + 1;
                    }
                    if actual < card {
                        bail!("tiny bitmap run cardinality exceeds its payload");
                    }
                    validate_tiny_runs(&data[..words], card)?;
                }
                _ => bail!("invalid tiny bitmap mode or cardinality"),
            }
            Ok(Bitmap {
                values: card as u64,
                maximum: None,
            })
        }
        _ => bail!("invalid Bifrost bitcontainer tag"),
    }
}

// Match the current engine's RLE cursor advances, including its look-ahead after
// the final emitted run. A superficially plausible but short payload can otherwise
// panic in that decoder. Inspect bounded words without allocating decoded values.
fn validate_tiny_runs(data: &[u16], card: usize) -> Result<()> {
    if card == 0 {
        return Ok(());
    }
    let offset = u32::from(data[2]) << 16;
    let (mut i, mut j, mut emitted) = (3, 4, 0);
    let mut value = offset | u32::from(*data.get(i).context("truncated tiny RLE start")?);
    while i < card + 3 && emitted < card {
        emitted += 1;
        value = value
            .checked_add(1)
            .context("tiny RLE cursor overflows the engine's u32 representation")?;
        if value & 0xffff > u32::from(*data.get(j).context("truncated tiny RLE end")?) {
            i += 2;
            j += 2;
            if i >= card + 3 {
                break;
            }
            value = offset | u32::from(*data.get(i).context("truncated tiny RLE look-ahead")?);
        }
    }
    Ok(())
}
