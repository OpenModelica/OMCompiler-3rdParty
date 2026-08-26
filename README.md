# OMCompiler-3rdParty

Third party sources used by OMCompiler.

## Known issue: SuiteSparse always shows up as dirty

Building SuiteSparse regenerates `SuiteSparse/SuiteSparse_config/SuiteSparse_config.h`
in the *source* tree, and that file is tracked by git. So after every build the
SuiteSparse submodule is dirty, which in turn makes this repo and the top-level
OpenModelica repo report the submodule as modified.

Typically the only difference is the timer macro, because we build
SuiteSparse_config without OpenMP:

```diff
-    #define SUITESPARSE_CONFIG_TIMER omp_get_wtime
+    #define SUITESPARSE_CONFIG_TIMER clock_gettime
```

The cause is `SuiteSparse_config/CMakeLists.txt`, which does a `configure_file`
into `${PROJECT_SOURCE_DIR}` instead of the build directory.

This is *not* going to be fixed upstream. See
[SuiteSparse#332](https://github.com/DrTimothyAldenDavis/SuiteSparse/issues/332)
and [SuiteSparse#340](https://github.com/DrTimothyAldenDavis/SuiteSparse/pull/340):
a configured copy of the header has to stay in the repository for the MATLAB
`mex` build path, which has no CMake or Makefile to generate one. Please don't
open a new ticket for it, and never commit the regenerated header.

### Workaround

Set the `skip-worktree` bit on that one file to hide the change from git:

```bash
git -C SuiteSparse update-index --skip-worktree SuiteSparse_config/SuiteSparse_config.h
```

`git status` is then clean in SuiteSparse, here, and at the OpenModelica top
level. The file on disk is untouched, the build still uses the regenerated
values.

Notes:

* The flag is local to your clone (it lives in `.git/modules/.../index`), so it
  cannot be committed for everyone and is lost when the submodule is re-cloned
  or `deinit`ed. CI is unaffected.
* List what is flagged with `git -C SuiteSparse ls-files -v | grep '^S'`.
* Undo with
  `git -C SuiteSparse update-index --no-skip-worktree SuiteSparse_config/SuiteSparse_config.h`.
* **Bumping the SuiteSparse submodule needs the flag removed first.** If the
  header changed in the new version, `git checkout` refuses with "local changes
  would be overwritten":

  ```bash
  git -C SuiteSparse update-index --no-skip-worktree SuiteSparse_config/SuiteSparse_config.h
  git -C SuiteSparse checkout -- SuiteSparse_config/SuiteSparse_config.h
  # ... bump the submodule, then set the flag again
  ```

A blunter alternative is `git config submodule.SuiteSparse.ignore dirty`, which
hides *any* worktree change in SuiteSparse (submodule commit changes are still
reported). It survives submodule bumps, but also hides real edits you make there.
