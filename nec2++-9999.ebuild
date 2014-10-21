# Copyright 1999-2013 Gentoo Foundation
# Distributed under the terms of the GNU General Public License v2
# $Header: $

EAPI=4

inherit eutils autotools flag-o-matic git-2 

DESCRIPTION="C++ version of the NEC2 antenna analysis program"
HOMEPAGE="http://elec.otago.ac.nz/w/index.php/Necpp"
EGIT_PROJECT='necpp'
EGIT_REPO_URI="https://github.com/tmolteno/necpp.git"

LICENSE="GPL-2"
SLOT="0"
KEYWORDS="amd64"
IUSE=""

DEPEND="sys-devel/libtool"

src_prepare() {
    eautoreconf
}

src_configure() {
    econf --without-lapack || die "Configuration failed"
}

src_compile() {
    replace-flags -O3 -O2
    replace-flags -Os -O2
    # emake has trouble with -j > 1
    emake -j1
}

src_install() {
    dobin src/nec2diff src/nec2++ || die "Binaries installation failed"
}
