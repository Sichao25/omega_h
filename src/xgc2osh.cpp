#include <Omega_h_file.hpp>
#include <Omega_h_cmdline.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();

  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh-in (xgc basename)");
  cmdline.add_arg<std::string>("mesh-out.osh");

  if (!cmdline.parse_final(comm, &argc, argv)) return -1;
  auto mesh_in = cmdline.get<std::string>("mesh-in (xgc basename)");
  auto mesh_out = cmdline.get<std::string>("mesh-out.osh");

  auto mesh = Omega_h::xgc::read(mesh_in, comm);
  Omega_h::binary::write(mesh_out, &mesh);

  return 0;
}
