# bigwigs-to-multivec
Convert multiple BigWig files to a single multivec file to be used in [HiGlass](http://higlass.io/)'s [horizontal multivec tracks](https://docs.higlass.io/track_types.html#horizontal-multivec).

## How to Use
Aggregate multiple bigwig files to a multivec file:
```

python convert.py input_files.txt [1000 [output_file.multires.mv5]]
```
where `input_files.txt` contains the paths of BigWig files, each line corresponding to each input file (Refer to the `example` folder) and 
`1000` correspond to the starting resolution (`default=1`).

Upload the multivec output into the [HiGlass server](https://github.com/higlass/higlass-server):
```
python manage.py ingest_tileset --filename my.multivec.file --filetype multivec /
       --datatype multivec --project-name "Experiment 1" --coordSystem hg19
```

Use the following track definition in HiGlass:
```json
{
    "type": "horizontal-multivec",
    "uid": "cistrome-track-1",
    "tilesetUid": "your-tile-set-uid",
    "server": "http://localhost:8001/api/v1",
    "options": {
        "labelPosition": "hidden",
        "labelColor": "black",
        "labelTextOpacity": 0.4,
        "valueScaling": "log",
        "trackBorderWidth": 0,
        "trackBorderColor": "black",
        "heatmapValueScaling": "log",
        "name": "your.multires.mv5",
        "labelLeftMargin": 0,
        "labelRightMargin": 0,
        "labelTopMargin": 0,
        "labelBottomMargin": 0,
        "labelShowResolution": true,
        "minHeight": 100,
        "colorbarPosition": "bottomLeft",
        "labelShowAssembly": true,
        "colorbarBackgroundColor": "#ffffff",
        "scaleStartPercent": "0.00000",
        "scaleEndPercent": "1.00000",
        "selectedRow": null
    },
    "width": 2530,
    "height": 458,
    "resolutions": [
        16384000,
        8192000,
        4096000,
        2048000,
        1024000,
        512000,
        256000,
        128000,
        64000,
        32000,
        16000,
        8000,
        4000,
        2000,
        1000
    ]
}
```