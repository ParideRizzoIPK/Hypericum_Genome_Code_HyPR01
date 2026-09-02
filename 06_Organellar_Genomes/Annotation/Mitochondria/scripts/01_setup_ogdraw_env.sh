#!/bin/bash
# ==============================================================================
# Reproducible setup for a local, offline OGDraw rendering engine.
#
# This builds a mamba environment (`ogdraw_env`) that can run the REAL
# OGDraw engine locally -- not an approximation. GeneMap-1.1.1
# (https://chlorobox.mpimp-golm.mpg.de/OGDraw-Downloads.html) is the actual
# Perl codebase behind the OGDraw web service (GeneMap/Linear.pm is the
# exact file that appeared in error logs from the web tool run earlier in
# this project -- same code, confirmed).
#
# Prerequisite: download GeneMap-1.1.1 from the URL above and note its path
# (referred to below as GENEMAP_DIR).
#
# Architecture note: this uses an x86_64 (osx-64) environment run under
# Rosetta 2 on Apple Silicon, NOT native arm64. Several required packages
# (perl-bioperl, perl-gd, the specific compiler wrapper Perl's build expects)
# only have osx-64 builds on bioconda/conda-forge -- there is no arm64 build
# available for this exact combination as of this writing. This is a
# deliberate, working choice, not an oversight.
# ==============================================================================
set -e -u -o pipefail

ENV_NAME="ogdraw_env"

echo "=== Step 1: BioPerl + GD (the two heaviest, version-sensitive deps) ==="
echo "    Must be created together in one solve -- installing separately risked"
echo "    a conflicting perl version being pulled in (this happened once during"
echo "    development: installing perl-app-cpanminus afterward silently pulled"
echo "    in perl 5.26.2 and broke the working 5.32.1 install; environment had"
echo "    to be removed and rebuilt). Recreate from scratch if in doubt, rather"
echo "    than adding packages to an existing solve."
CONDA_SUBDIR=osx-64 mamba create -n "$ENV_NAME" -y --channel-priority flexible \
  -c bioconda -c conda-forge perl-bioperl perl-gd

echo "=== Step 2: PostScript::Simple (conda-packaged, no conflict) ==="
CONDA_SUBDIR=osx-64 mamba install -n "$ENV_NAME" -y --channel-priority flexible \
  -c bioconda -c conda-forge perl-postscript-simple

echo "=== Step 3: ImageMagick C library + compiler toolchain (for Image::Magick) ==="
CONDA_SUBDIR=osx-64 mamba install -n "$ENV_NAME" -y --channel-priority flexible \
  -c conda-forge -c bioconda imagemagick clang_osx-64 clangxx_osx-64

echo "=== Step 4: Bio::Perl.pm -- missing from this bioperl conda package (thin"
echo "    convenience wrapper, apparently trimmed from the bioconda build)."
echo "    Fetched directly from the closest matching bioperl-live release tag"
echo "    (release-1-7-2; the installed conda package is bioperl 1.7.8, which"
echo "    has no matching git tag -- Bio::Perl.pm is stable across these minor"
echo "    versions so this is safe)."
ENV_PREFIX=$(mamba run -n "$ENV_NAME" perl -e 'print $^X' | sed 's|/bin/perl||')
BIO_DIR=$(mamba run -n "$ENV_NAME" perl -MBio::SeqIO -e 'print $INC{"Bio/SeqIO.pm"}' | sed 's|SeqIO.pm||')
curl -sL "https://raw.githubusercontent.com/bioperl/bioperl-live/release-1-7-2/Bio/Perl.pm" \
  -o "${BIO_DIR}Perl.pm"

echo "=== Step 5: Bio::Restriction::* -- also missing from the conda package."
echo "    Fetched the whole subtree (Analysis.pm + its own nested deps) from"
echo "    the same release tag."
curl -sL "https://github.com/bioperl/bioperl-live/archive/refs/tags/release-1-7-2.tar.gz" \
  -o /tmp/bioperl-full.tar.gz
tar xzf /tmp/bioperl-full.tar.gz -C /tmp
cp -r /tmp/bioperl-live-release-1-7-2/Bio/Restriction "${BIO_DIR}../Bio/"

echo "=== Step 6: XML::Generator, XML::Simple -- pure-Perl, not on bioconda."
echo "    Installed via cpanm, but invoked using the conda env's OWN perl"
echo "    interpreter explicitly (not a conda-installed cpanm -- that's what"
echo "    caused the Step 1 breakage) so the modules land in the right place"
echo "    without touching the environment's package solve."
"${ENV_PREFIX}/bin/perl" /opt/homebrew/bin/cpanm --notest XML::Generator XML::Simple

echo "=== Step 7: Image::Magick XS bindings -- built against the conda-env's"
echo "    own ImageMagick library (not any system/Homebrew one)."
export PKG_CONFIG_PATH="${ENV_PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export PATH="${ENV_PREFIX}/bin:$PATH"
export CPATH="${ENV_PREFIX}/include:${ENV_PREFIX}/include/ImageMagick-7"
export LIBRARY_PATH="${ENV_PREFIX}/lib"
"${ENV_PREFIX}/bin/perl" /opt/homebrew/bin/cpanm --notest Image::Magick

echo
echo "=== Verifying full dependency chain ==="
mamba run -n "$ENV_NAME" perl -v | head -3
for mod in Bio::Perl Bio::Restriction::Analysis Bio::SeqIO GD Image::Magick \
           PostScript::Simple XML::Generator XML::Simple; do
    mamba run -n "$ENV_NAME" perl -M"$mod" -e "print '$mod OK\n';" 2>&1
done

echo
echo "Setup complete. Test with:"
echo "  cd <GENEMAP_DIR> && mamba run -n $ENV_NAME perl -Ilib bin/drawgenemap --help"
