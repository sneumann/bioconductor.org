Mirrors
=======================================

If you are interested in maintaining a mirror of this site (for either
public or private use) [read this](mirror-how-to/). You can check current mirror
status on the [Bioconductor dashboard](/dashboard/) mirror section. To select a
Bioconductor mirror use the R function `chooseBioCmirror()`

<% for country in config[:mirrors] %>
  <% next if country.keys.first.to_s == "0-Bioconductor" %>
<%= country.keys.first.to_s %>
------------------------

<% for mirror in country.values.first %>
* [<%= mirror[:institution] %>](<%= mirror[:institution_url] %>)

  URLs: <%= render_mirror_urls(mirror) %>

  Contact: <%= render_mirror_contacts(mirror) %>

<% end %>
<% end %>
