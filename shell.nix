{ pkgs ? import <nixpkgs> {}}:
pkgs.mkShell {
  packages = with pkgs; [
    libgcc
    onetbb
    hwloc
    cmake
    gdb
  ];

  buildInputs = [ pkgs.bashInteractive ];
}