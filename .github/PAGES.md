# Package website

The `pkgdown` workflow builds the checked-out package and its documentation on
Ubuntu. Pull requests targeting `master` build the site and upload a reviewable
`github-pages` artifact. Pushes to `master` also deploy that artifact to
<https://gestimation.github.io/cifmodeling/> after the build succeeds.

## One-time migration

1. Confirm that the workflow's **Build pkgdown site** check passes on the PR.
2. In **Settings > Pages > Build and deployment > Source**, select
   **GitHub Actions** instead of **Deploy from a branch**. Keep the existing
   site URL and domain settings.
3. Merge the workflow PR into `master`. Its push run builds and deploys the
   current site. If it was merged before changing the Pages setting, open
   **Actions > pkgdown > Run workflow**, select `master`, and run it again.
4. Confirm the **Deploy pkgdown site** job succeeds, then check the public
   version number, Reference pages and both reporting articles.

The checked-in `docs/` directory is the old publication snapshot. The new
workflow generates its complete site in a fresh runner-temporary directory;
it neither publishes the snapshot nor commits generated files. Leave `docs/`
in place during this migration; removing historical generated files can be a
separate cleanup after successful publication.
Historical article URLs are retained through redirects in `_pkgdown.yml`.

## Routine updates

Update R documentation, `vignettes/`, NEWS, or `_pkgdown.yml` in a normal PR.
Merge into `master` after review; no local site build or HTML commit is needed.
The root site follows `master` (release display mode), including when the
package version has not changed. Manual runs on other branches build only.

The build installs the local package and its suggested/website dependencies,
runs the examples and articles, and checks that every vignette has an article
page and that both gtsummary articles contain rendered tables. Only the deploy
job has Pages write and OIDC permissions; PRs do not deploy. Production runs
are serialized, and a failed build leaves the last successful site online.

For a failed build, inspect the pkgdown job log and fix the source in a PR.
For a failed deployment, first check the Pages source above and the
`github-pages` environment's deployment restrictions. The workflow uses the
standard `GITHUB_TOKEN`; no personal access token or deploy key is needed.

References: [GitHub Pages custom workflows](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages),
[r-lib pkgdown workflow](https://github.com/r-lib/actions/blob/v2/examples/pkgdown.yaml).
