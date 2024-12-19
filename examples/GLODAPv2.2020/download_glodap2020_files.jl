import Downloads
import ZipFile

function download_glodap2020_files(;
    urlroot="https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0210813/",
    matfiles=[
        "GLODAPv2.2020_Arctic_Ocean.mat",
        "GLODAPv2.2020_Atlantic_Ocean.mat",
        "GLODAPv2.2020_Indian_Ocean.mat",
        "GLODAPv2.2020_Pacific_Ocean.mat",
    ]
)

    for matfile in matfiles
        url = urlroot*matfile
        @info "Downloading $url -> $matfile"
        Downloads.download(url, matfile)
    end

    return nothing
end