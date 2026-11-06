process FETCH_MITOS_REFDB {
    tag 'mitos_refdb'
    label 'process_single'

    // Plain Python download: Nextflow's own HTTP staging fails on Zenodo URLs on some
    // hosts ("Unable to access path"), while a normal client works.
    conda "${moduleDir}/../circos_plot/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pycirclize_python:731acb4cc6a99f3c' :
        'community.wave.seqera.io/library/pycirclize_python:ad50b4c46d2a3f7e' }"

    input:
    val url

    output:
    path "mitos_refdb.tar.bz2", emit: archive

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python - <<'PY'
    import shutil, time, urllib.request
    for attempt in range(5):
        try:
            with urllib.request.urlopen("${url}", timeout=300) as r, open("mitos_refdb.tar.bz2", "wb") as fh:
                shutil.copyfileobj(r, fh)
            break
        except Exception as err:
            print(f"download attempt {attempt + 1} failed: {err}")
            time.sleep(10 * (attempt + 1))
    else:
        raise SystemExit("could not download ${url}")
    PY
    """

    stub:
    """
    touch mitos_refdb.tar.bz2
    """
}
