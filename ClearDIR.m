function ClearDIR
   GROPTBase=pwd;
   fprintf('Deleting GROPT paths from %s...\n',GROPTBase);
   rmpath(genpath(GROPTBase));
end