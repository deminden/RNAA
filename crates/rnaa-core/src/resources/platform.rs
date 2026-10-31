#[cfg(unix)]
use anyhow::Context;
use anyhow::Result;
use std::path::Path;

pub(super) fn memory_capacity_bytes() -> Result<u64> {
    #[cfg(target_os = "macos")]
    {
        let mut memory = 0_u64;
        let mut size = std::mem::size_of::<u64>();
        // SAFETY: the name is NUL-terminated; output is writable u64 storage
        // with its exact size. No new sysctl value is supplied.
        if unsafe {
            libc::sysctlbyname(
                c"hw.memsize".as_ptr(),
                (&mut memory as *mut u64).cast(),
                &mut size,
                std::ptr::null_mut(),
                0,
            )
        } != 0
        {
            return Err(std::io::Error::last_os_error())
                .context("cannot detect host memory capacity");
        }
        if size != std::mem::size_of::<u64>() || memory == 0 {
            anyhow::bail!("invalid host memory capacity; configure resources.memory_mb");
        }
        Ok(memory)
    }
    #[cfg(target_os = "linux")]
    {
        // SAFETY: sysconf takes a constant selector and no pointer arguments.
        let pages = unsafe { libc::sysconf(libc::_SC_PHYS_PAGES) };
        let page_size = unsafe { libc::sysconf(libc::_SC_PAGESIZE) };
        let mut memory = u64::try_from(pages)
            .ok()
            .zip(u64::try_from(page_size).ok())
            .and_then(|(a, b)| a.checked_mul(b))
            .filter(|v| *v > 0)
            .context("cannot detect host memory; configure resources.memory_mb")?;
        for path in [
            "/sys/fs/cgroup/memory.max",
            "/sys/fs/cgroup/memory/memory.limit_in_bytes",
        ] {
            if let Ok(value) = std::fs::read_to_string(path)
                && let Ok(limit) = value.trim().parse::<u64>()
                && limit > 0
            {
                memory = memory.min(limit);
            }
        }
        Ok(memory)
    }
    #[cfg(not(any(target_os = "macos", target_os = "linux")))]
    anyhow::bail!("automatic memory detection requires macOS/Linux; configure resources.memory_mb")
}

#[cfg(unix)]
pub(super) fn available_disk_bytes(path: &Path) -> Result<u64> {
    use std::{ffi::CString, mem::MaybeUninit, os::unix::ffi::OsStrExt};
    let path =
        CString::new(path.as_os_str().as_bytes()).context("filesystem path contains a NUL byte")?;
    let mut info = MaybeUninit::<libc::statvfs>::uninit();
    // SAFETY: path is NUL-terminated and info points to writable, correctly sized
    // storage. A successful statvfs call initializes the entire output structure.
    if unsafe { libc::statvfs(path.as_ptr(), info.as_mut_ptr()) } != 0 {
        return Err(std::io::Error::last_os_error()).context("cannot inspect available disk space");
    }
    // SAFETY: statvfs succeeded above.
    let info = unsafe { info.assume_init() };
    let block = if info.f_frsize == 0 {
        info.f_bsize
    } else {
        info.f_frsize
    };
    let bytes = u128::from(info.f_bavail) * u128::from(block);
    u64::try_from(bytes).context("available disk space exceeds supported byte accounting")
}

#[cfg(not(unix))]
pub(super) fn available_disk_bytes(_path: &Path) -> Result<u64> {
    anyhow::bail!("disk reservation probing is supported on Unix hosts")
}
