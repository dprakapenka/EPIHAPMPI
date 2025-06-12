-- rewrite.links.lua

-- Cache mapping from file path to heading ID
local heading_map = {}

-- turn arbitrary heading text into the same ID Pandoc generates:
local function slugify(text)
  text = text:lower()
  -- strip punctuation except hyphens/spaces
  text = text:gsub("[^%w%s%-.]", "")
  -- spaces to hyphens
  text = text:gsub("%s+", "-")
  return text
end

function Link(el)
  if el.target and el.target:match("[^/]+%.md$") then
    local path = el.target
    io.stderr:write(("DEBUG: original target='%s'\n"):format(path))

    local id = heading_map[path]
    if not id then
      local f = io.open(path, "r")
      if f then
        for line in f:lines() do
          local text = line:match("^%s*#+%s*(.+)")
          if text then
            id = slugify(text)
            heading_map[path] = id
            break
          end
        end
        f:close()
      end
      if not id then
        -- fallback to filename
        local name = path:match("([^/]+)%.md$")
        id = slugify(name)
        heading_map[path] = id
      end
    end

    io.stderr:write(("DEBUG: resolved slugified id='%s'\n"):format(id))
    el.target = "#" .. id
    io.stderr:write(("DEBUG: new target='%s'\n"):format(el.target))
  end
  return el
end
