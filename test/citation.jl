import Test
import TOML

function test_citation_metadata()
    Test.@testset "Citation metadata" begin
        repo_root = normpath(joinpath(@__DIR__, ".."))
        read_text(path) = replace(read(path, String), "\r\n" => "\n")

        project = TOML.parsefile(joinpath(repo_root, "Project.toml"))
        citation = read_text(joinpath(repo_root, "CITATION.cff"))
        readme = read_text(joinpath(repo_root, "README.md"))
        checklist =
            read_text(joinpath(repo_root, ".github", "RELEASE_CHECKLIST.md"))

        concept_doi = "10.5281/zenodo.7511329"
        version_doi = "10.5281/zenodo.20926960"
        upstream_doi = "10.1287/ijoc.2017.0798"
        ecosystem_doi = "10.1080/10556788.2026.2702926"
        version = string(project["version"])

        Test.@test occursin("doi: \"$concept_doi\"", citation)
        Test.@test occursin("value: \"$version_doi\"", citation)
        Test.@test occursin("version: \"$version\"", citation)
        Test.@test occursin(
            "description: \"Zenodo DOI for version $version\"",
            citation,
        )
        Test.@test occursin("upstream MQLib", citation)

        Test.@test occursin("badge/DOI/$concept_doi.svg", readme)
        Test.@test occursin("doi.org/$concept_doi", readme)
        Test.@test occursin("doi.org/$version_doi", readme)
        Test.@test occursin("v$version", readme)
        Test.@test occursin(upstream_doi, readme)
        Test.@test occursin(ecosystem_doi, readme)
        Test.@test occursin(
            "does not identify or archive upstream MQLib",
            readme,
        )

        Test.@test occursin(concept_doi, checklist)
        Test.@test occursin("Can manage", checklist)
        Test.@test occursin(
            "cffconvert --validate --infile CITATION.cff",
            checklist,
        )
    end

    return nothing
end
