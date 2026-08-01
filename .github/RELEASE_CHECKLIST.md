# Release checklist

Use this checklist for every tagged release. The Zenodo concept DOI
`10.5281/zenodo.7511329` is the permanent identifier for the MQLib.jl wrapper;
do not create another concept record for a new version. It does not identify
or archive the upstream [MQLib](https://github.com/MQLib/MQLib) project.

## Before tagging

- [ ] Confirm at least two active JuliaQUBO maintainers have **Can manage**
      access to the Zenodo record and all future versions. Use separate personal
      accounts; do not share a Zenodo account, password, or API token.
- [ ] Update the version in `Project.toml`, the release date and version in
      `CITATION.cff`, and the release notes in `CHANGELOG.md`.
- [ ] Keep the concept DOI in `CITATION.cff` and the README badge unchanged.
- [ ] Confirm the wrapper creators and their order, affiliations, and ORCIDs.
- [ ] Validate the citation metadata with
      `cffconvert --validate --infile CITATION.cff`.
- [ ] Run the package test suite with
      `julia --project -e 'using Pkg; Pkg.test()'`.

## After publishing the GitHub release

- [ ] Confirm the Zenodo GitHub integration published a **new version** under
      the existing concept record automatically when TagBot created the GitHub
      release. Do not deposit a new upload or concept DOI by hand; if the
      integration is ever disabled, re-enable it rather than archiving
      manually.
- [ ] Confirm the archived Zenodo record matches the official GitHub release's
      tag and version, MIT license, repository URL, Julia package UUID
      (`16f11440-1623-44c9-850c-358a6c72f3c9`), creators, affiliations, and
      ORCIDs.
- [ ] Confirm the description and related identifiers attribute upstream MQLib
      without implying that the wrapper's DOI archives upstream MQLib.
- [ ] Record the published Zenodo version DOI in the GitHub release notes.
      Keep the concept DOI in `CITATION.cff`; Zenodo assigns the version DOI
      only after it processes the GitHub release, so do not pin a version DOI
      in the release archive's CFF.
- [ ] Verify the version DOI resolves to that exact archive and the concept DOI
      resolves to the latest archived release:

      ```sh
      curl --fail --location --output /dev/null https://doi.org/<version-doi>
      curl --fail --location --output /dev/null https://doi.org/10.5281/zenodo.7511329
      curl --fail --location https://zenodo.org/api/records/<version-record-id>
      ```

      In the API response, confirm `conceptdoi` is
      `10.5281/zenodo.7511329`, `metadata.version` matches the release tag,
      `metadata.custom["code:codeRepository"]` identifies this repository,
      and the record is the latest version under the concept. If Zenodo responds
      slowly or returns HTTP 429, wait for its `Retry-After` interval and retry
      once before treating the DOI as broken.

- [ ] Download the Zenodo archive, record and verify its checksum, confirm it
      contains the tagged `Project.toml` version, and verify its source commit
      matches the GitHub release tag.
