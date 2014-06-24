TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += src apps

OTHER_FILES += .gitignore README.md LICENSE.txt


DIRS = scratch/QMC_SCRATCH/walker_positions

for(DIR, DIRS) {
     mkcommands += $$shadowed($$PWD)/$$DIR
}

createDirs.commands = $(MKDIR) $$mkcommands

first.depends = $(first) createDirs
export(first.depends)
export(createDirs.commands)

QMAKE_EXTRA_TARGETS += first createDirs

