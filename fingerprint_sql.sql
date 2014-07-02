select chromosome, position, length, count(distinct(length)) as n from variants where chromosome in ('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8') and length != 0 and Oxford not like '%N%' and Oxford not like '%Y%' and Oxford not like '%M%' and Oxford not like '%K%' and Oxford not like '%R%' and Oxford not like '%W%' and Oxford not like '%S%' and consensus not like '%N%' and consensus not like '%Y%' and consensus not like '%M%' and consensus not like '%K%' and consensus not like '%R%' and consensus not like '%W%' and consensus not like '%S%' and abs(length) > 19 group by position limit 10;

select chromosome, position, length, count(distinct(length)) as n into outfile '/tmp/ma_indel.csv' fields terminated by ',' from variants where chromosome in ('Chr1','Chr2','Chr3','Chr4','Chr5','Chr6','Chr7','Chr8') and length != 0 and Oxford not like '%N%' and Oxford not like '%Y%' and Oxford not like '%M%' and Oxford not like '%K%' and Oxford not like '%R%' and Oxford not like '%W%' and Oxford not like '%S%' and consensus not like '%N%' and consensus not like '%Y%' and consensus not like '%M%' and consensus not like '%K%' and consensus not like '%R%' and consensus not like '%W%' and consensus not like '%S%' and abs(length) > 19 group by position order by n desc;