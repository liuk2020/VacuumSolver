pip uninstall vacuum
rm -rf build
meson setup build
meson compile -C build
pip install .
