# first running to load DLL if windows

.onLoad <-function(lib, pkg) {
    if (.Platform$OS.type == "windows") {
        arch_type <- .Platform$r_arch
        print(arch_type)
        dll_path <- system.file("libs/x64/rgeoda.dll", package="rgeoda")
        if (arch_type == "i386") {
            dll_path <- system.file("libs/i386/rgeoda.dll", package="rgeoda")
        } 
        dyn.load(dll_path)
        cacheMetaData(1)
    } else {
        # unix
        dll_path <- system.file("libs/rgeoda.so", package="rgeoda")
        dyn.load(dll_path)
        cacheMetaData(1)
    }
}