ZAF for macOS (Apple Silicon)
================================

Starting ZAF
------------
Double-click ZAF.app. If macOS blocks the first launch because the app is not
notarized, right-click ZAF.app, choose Open, and confirm that you want to open
it.

TEM instrument settings
-----------------------
ZAF creates the active instrument settings file when the app starts for the
first time. Its location is:

~/Library/Application Support/ZAF/ZAF_instrument_settings.txt

To open that folder in Finder:

1. Launch ZAF at least once so the file is created.
2. Open Finder.
3. Choose Go > Go to Folder from the menu bar at the top of the screen.
   You can also press Shift-Command-G.
4. Enter this path and press Return:

   ~/Library/Application Support/ZAF

Edit only the documented numeric values in ZAF_instrument_settings.txt, save
the file, and restart ZAF. ZAF preserves this per-user file when the app is
moved or upgraded.

Why the file is not beside ZAF.app
----------------------------------
macOS Gatekeeper may run a downloaded application from a temporary
App Translocation path. A settings file beside the original app may then be
invisible to the running program. The Application Support location avoids
that problem and prevents edits from being lost during upgrades.

