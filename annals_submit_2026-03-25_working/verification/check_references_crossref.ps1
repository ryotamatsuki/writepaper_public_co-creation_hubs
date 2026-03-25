$refs = @(
  @{Key='AgrawalCockburn2003'; Doi='10.1016/S0167-7187(03)00081-X'},
  @{Key='AsheimLawtonSmithOughton2011'; Doi='10.1080/00343404.2011.596701'},
  @{Key='BatheltMalmbergMaskell2004'; Doi='10.1191/0309132504ph469oa'},
  @{Key='Boschma2005'; Doi='10.1080/0034340052000320887'},
  @{Key='ChaminadeBellandiPlecheroSantini2019'; Doi='10.1080/09654313.2019.1610727'},
  @{Key='CoenenAsheimBuggeHerstad2017'; Doi='10.1177/0263774X16646583'},
  @{Key='CohenLevinthal1990'; Doi='10.2307/2393553'},
  @{Key='CookeUrangaEtxebarria1997'; Doi='10.1016/S0048-7333(97)00025-5'},
  @{Key='DurantonPuga2004'; Doi='10.1016/S1574-0080(04)80005-1'},
  @{Key='JaffeTrajtenbergHenderson1993'; Doi='10.2307/2118401'},
  @{Key='KlineMoretti2014ARE'; Doi='10.1146/annurev-economics-080213-041024'},
  @{Key='KlineMoretti2014QJE'; Doi='10.1093/qje/qjt034'},
  @{Key='MerrellPhillipsonGortonCowie2022'; Doi='10.1016/j.jrurstud.2022.05.016'},
  @{Key='MerrellRoweCowieGkartzios2021'; Doi='10.1177/02690942221085498'},
  @{Key='NeumarkSimpson2015'; Doi='10.1016/B978-0-444-59531-7.00018-1'},
  @{Key='OECD2025'; Doi='10.1787/c86de0f4-en'},
  @{Key='RocheOettlCatalini2024'; Doi='10.1287/mnsc.2022.03555'},
  @{Key='RosenthalStrange2004'; Doi='10.1016/S1574-0080(04)80006-3'},
  @{Key='StorperVenables2004'; Doi='10.1093/jeg/4.4.351'},
  @{Key='TodtlingTrippl2005'; Doi='10.1016/j.respol.2005.01.018'},
  @{Key='TorreRallet2005'; Doi='10.1080/0034340052000320842'},
  @{Key='VoglAkhavan2022'; Doi='10.1108/JPIF-12-2021-0108'}
)

$ProgressPreference = 'SilentlyContinue'

$results = foreach ($ref in $refs) {
  $uri = 'https://api.crossref.org/works/' + [uri]::EscapeDataString($ref.Doi)
  try {
    $msg = (Invoke-RestMethod -UseBasicParsing -Uri $uri).message
    $authors = @(
      $msg.author | ForEach-Object {
        if ($_.family) {
          "$($_.family), $($_.given)"
        } elseif ($_.name) {
          $_.name
        }
      }
    ) -join '; '

    [PSCustomObject]@{
      Key = $ref.Key
      DOI = $ref.Doi
      Title = (@($msg.title) -join ' | ')
      Container = (@($msg.'container-title') -join ' | ')
      Year = if ($msg.issued.'date-parts') { $msg.issued.'date-parts'[0][0] } else { $null }
      Volume = $msg.volume
      Issue = $msg.issue
      Page = $msg.page
      Authors = $authors
    }
  } catch {
    [PSCustomObject]@{
      Key = $ref.Key
      DOI = $ref.Doi
      Title = 'ERROR'
      Container = $_.Exception.Message
      Year = $null
      Volume = $null
      Issue = $null
      Page = $null
      Authors = $null
    }
  }
}

$results | ConvertTo-Json -Depth 4
