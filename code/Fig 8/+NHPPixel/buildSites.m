function sites = buildSites(mgr)
    % spreadsheet in https://docs.google.com/spreadsheets/d/1s-IWQ3TMdFjl5mGkyC44DAvi5FHIz2mEYPPwBc9Aqbo/edit#gid=0

    n = 1;
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, subject="Cousteau", date=datetime("2021-03-18"), save_tags = 0, ...
        gNum=0, tNum=0, imecNum=0, npix_region="PMd");

    n = n + 1;
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, subject="Cousteau", date=datetime("2021-03-29"), save_tags=[0 1 2], ...
        sync_chan = 2, gNum=0, tNum=0, imecNum=0, npix_region="PMd");

    n = n + 1; 
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, subject="Cousteau", date=datetime("2021-05-07"), save_tags=0, ...
        sync_chan = 2, gNum=0, tNum=0, imecNum=0, npix_region="PMd");

    n = n + 1; 
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, subject="Cousteau", date=datetime("2021-05-20"), save_tags=0, ...
        sync_chan = 2, gNum=0, tNum=0, imecNum=0, npix_region="PMd");

    n = n + 1; 
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, subject="Igor", date=datetime("2022-05-18"), save_tags=1:3, ...
        lf_rmsRange=[1 80], gNum=0, tNum=0, imecNum=0, npix_region="M1");
    
    n = n + 1; 
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, subject="Igor", date=datetime("2022-05-19"), save_tags=3, ...
        gNum=0, tNum=0, imecNum=0, npix_region="M1");
    
    n = n + 1; 
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, subject="Igor", date=datetime("2022-05-20"), save_tags=1:5, ...
        gNum=0, tNum=0, imecNum=0, npix_region="M1");

    n = n + 1; 
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, subject="Igor", date=datetime("2022-05-25"), save_tags=0:2, ...
        gNum=0, tNum=0, imecNum=1, npix_region="M1");

    n = n + 1; 
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, id="Igor_20220506a", subject="Igor", date=datetime("2022-05-26"), save_tags=1, ...
        lf_rmsRange=[50 250], gNum=0, tNum=0, imecNum=0, npix_region="M1");

    n = n + 1; 
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, id="Igor_20220506b", subject="Igor", date=datetime("2022-05-26"), save_tags=1, ...
        gNum=0, tNum=0, imecNum=1, npix_region="M1");

    n = n + 1; 
    sites(n) = NHPPixel.PacmanSessionInfo(mgr, id="Igor_20220506c", subject="Igor", date=datetime("2022-05-26"), save_tags=1, ...
        gNum=0, tNum=0, imecNum=2, npix_region="M1");
end